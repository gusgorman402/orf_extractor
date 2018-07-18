#!/usr/bin/perl -w

use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::SeqIO;
use strict;

my $report = Bio::SearchIO->new( -format => 'blast', -file => $ARGV[0] );
my $cdna_file = Bio::DB::Fasta->new( $ARGV[1] );
my $prot_file = Bio::SeqIO->new( -file => ">".$ARGV[1].".prot", -format => 'fasta' );

my $hit_counter;
my %one_hsp_transcript = ();

while( my $result = $report->next_result )
{
    if( $result->num_hits > 0 )
    {
        $hit_counter++;
        my $hit = $result->next_hit;
        if( $hit->num_hsps == 1 )
        {
            my $hsp = $hit->next_hsp;
            my $best_hit_id = $hit->name;
            if( !defined( $one_hsp_transcript{$best_hit_id} ) 
                || $hsp->score > $one_hsp_transcript{$best_hit_id} )
            {
                $one_hsp_transcript{$best_hit_id} = $hsp->score;
            }

            my $blast_orf = $hsp->query_string;
            $blast_orf =~ s/\-//g;
            if( $blast_orf =~ /^(\w+)\*.+/ )
            {
                if( ($hsp->end('hit') == $hit->length) && (length($blast_orf) - length($1) < 10) )
                {
                    $blast_orf = $1;
                }
                else
                {
                    next;
                }
            }

            my $query_id = $result->query_name;
            my $cdna_seq = $cdna_file->get_Seq_by_id($query_id);
            my ($long_orf, $orf_start, $orf_stop );

            if( $cdna_seq->seq =~ /N/ )
            {
                $long_orf = $blast_orf;
                if( $hsp->strand('query') < 1 )
                {
                    $orf_start = $hsp->end('query');
                    $orf_stop = $hsp->start('query');
                }
                else
                {
                    $orf_start = $hsp->start('query');
                    $orf_stop = $hsp->end('query');
                }
            }
            else
            {
                my $temp_input_file = Bio::SeqIO->new( -file => ">temp".$hit_counter.".fna", -format => 'fasta' );
                $temp_input_file->write_seq( $cdna_seq );
            
                my $cmd = "getorf -sequence temp".$hit_counter.".fna -outseq temp".$hit_counter.".orf -minsize 90 -auto";
            
                system( $cmd );
                my $temp_orf_file = Bio::SeqIO->new( -file => "temp".$hit_counter.".orf", -format => 'fasta' );

                while( my $temp_seq = $temp_orf_file->next_seq )
                {
                    if( $temp_seq->seq =~ m/$blast_orf/ )
                    {
                        $long_orf = $temp_seq->seq;
                        $temp_seq->desc =~ /.*\[(\d+)\s\-\s(\d+)\].+/;
                        $orf_start = $1;
                        $orf_stop = $2;
                        last;
                    }
                }
                $cmd = "rm temp".$hit_counter.".*";
                system( $cmd );
                if( $long_orf !~ /\w+/ ){ next }
            }
            
            my $start_codon = "NO";
            if( $hsp->start('hit') == 1 && $hsp->hit_string =~ /^M.+/ )
            {
                if( $blast_orf =~ /^M.+/ )
                {
                    $start_codon = "YES";
                    if( $long_orf =~ /(\w+)($blast_orf.*)/ )
                    {
                        my $extra_dna = length( $1 ) * 3;
                        $long_orf = $2;
                        if( $orf_start > $orf_stop ){ $orf_start = $orf_start - $extra_dna }
                        elsif( $orf_start != 0 ){ $orf_start = $orf_start + $extra_dna }
                    }
                }
            }
            elsif( $hsp->start('hit') < 30 )
            {
                if( $long_orf =~ /^(\w*)(M\w*$blast_orf.*)/ )
                {
                    $start_codon = "YES";
                    my $extra_dna = length( $1 ) * 3;
                    $long_orf = $2; 
                    if( $orf_start > $orf_stop ){ $orf_start = $orf_start - $extra_dna }
                    elsif( $orf_start != 0 ){ $orf_start = $orf_start + $extra_dna }
                }
            }

            my $stop_codon = "YES";
            if( $orf_start < $orf_stop && $orf_start < 5 && $cdna_seq->length - $orf_stop < 5 ){ $stop_codon = "NO" }
            if( $orf_start > $orf_stop && $orf_stop < 4 && $cdna_seq->length - $orf_start < 5 ){ $stop_codon = "NO" }
            if( $stop_codon eq "YES" && $hsp->end('hit') / $hit->length > 0.85 )
            {
                if( $orf_start < $orf_stop && abs( $orf_stop - $hsp->end('query')) < 50){ $stop_codon = "YES_REAL" }
                if( $orf_start > $orf_stop && abs( $orf_start - $hsp->end('query')) < 50){ $stop_codon = "YES_REAL" }
            }

            my $new_desc = length($long_orf).":".$orf_start."-".$orf_stop.":".length($blast_orf).":".$hsp->start('query')."-".$hsp->end('query');
            $new_desc = $new_desc.":".$hit->name.":".$hit->length.":".$hsp->start('hit')."-".$hsp->end('hit').":".$start_codon.":".$stop_codon;

            my $new_seq = Bio::Seq->new( -display_id => $query_id,
                                         -desc => $new_desc,
                                         -seq => $long_orf );
            $prot_file->write_seq( $new_seq );
        }
    }
}

my $repeat_report = Bio::SearchIO->new( -format => 'blast', -file => $ARGV[0] );

while( my $result = $repeat_report->next_result )
{
    my $long_orf = "";
    if( $result->num_hits > 0 )
    {
        my $hit = $result->next_hit;
        if( $hit->num_hsps > 1 )
        {
            $hit_counter++;
            my $hsp = $hit->next_hsp;
            my $hsp_status = "MultiHSP";
            if( !defined( $one_hsp_transcript{$hit->name} ) || $hsp->score > $one_hsp_transcript{$hit->name} )
            {
                if( !defined( $one_hsp_transcript{$hit->name} )) { $hsp_status = "UniqueMultiHSP" }
                elsif( $hsp->score > $one_hsp_transcript{$hit->name} ){ $hsp_status = "BetterMultiHSP" }

                my $blast_orf = $hsp->query_string;
                $blast_orf =~ s/\-//g;
                if( $blast_orf =~ /^(\w+)\*.+/ )
                {
                    if( ($hsp->end('hit') == $hit->length) && (length($blast_orf) - length($1) < 10) )
                    {
                        $blast_orf = $1;
                    }
                    else
                    {
                        next;
                    }
                }

                my $query_id = $result->query_name;
                my $cdna_seq = $cdna_file->get_Seq_by_id( $query_id );
                my ($long_orf, $orf_start, $orf_stop );
                
                if( $cdna_seq->seq =~ /N/ )
                {
                    $long_orf = $blast_orf;
                    $orf_start = 0;
                    $orf_stop = 0;
                }
                else
                {
                    my $temp_input_file = Bio::SeqIO->new( -file => ">temp".$hit_counter.".fna", -format => 'fasta' );
                    $temp_input_file->write_seq( $cdna_seq );
            
                    my $cmd = "getorf -sequence temp".$hit_counter.".fna -outseq temp".$hit_counter.".orf -minsize 90 -auto";
            
                    system( $cmd );
                    my $temp_orf_file = Bio::SeqIO->new( -file => "temp".$hit_counter.".orf", -format => 'fasta' );

                    while( my $temp_seq = $temp_orf_file->next_seq )
                    {
                        if( $temp_seq->seq =~ m/$blast_orf/ )
                        {
                            $long_orf = $temp_seq->seq;
                            $temp_seq->desc =~ /.*\[(\d+)\s\-\s(\d+)\].+/;
                            $orf_start = $1;
                            $orf_stop = $2;
                            last;
                        }
                    }
                    $cmd = "rm temp".$hit_counter.".*";
                    system( $cmd );
                    if( $long_orf !~ /\w+/ ){ next }
                }

                my $start_codon = "NO";
                if( $hsp->start('hit') == 1 && $hsp->hit_string =~ /^M.+/ )
                {
                    if( $blast_orf =~ /^M.+/ )
                    {
                        $start_codon = "YES";
                        if( $long_orf =~ /(\w+)($blast_orf.*)/ )
                        {
                            my $extra_dna = length( $1 ) * 3;
                            $long_orf = $2;
                            if( $orf_start > $orf_stop ){ $orf_start = $orf_start - $extra_dna }
                            elsif( $orf_start != 0 ){ $orf_start = $orf_start + $extra_dna }
                        }
                    }
                }
                elsif( $hsp->start('hit') < 30 && $hsp->start('hit') > 1 )
                {
                    if( $long_orf =~ /^(\w*)(M\w*$blast_orf.*)/ )
                    {
                        $start_codon = "YES";
                        my $extra_dna = length( $1 ) * 3;
                        $long_orf = $2; 
                        if( $orf_start > $orf_stop ){ $orf_start = $orf_start - $extra_dna }
                        elsif( $orf_start != 0 ){ $orf_start = $orf_start + $extra_dna }
                    }
                }

                my $stop_codon = "YES";
                if( $orf_start < $orf_stop && $orf_start < 5 && $cdna_seq->length - $orf_stop < 5 ){ $stop_codon = "NO" }
                if( $orf_start > $orf_stop && $orf_stop < 4 && $cdna_seq->length - $orf_start < 5 ){ $stop_codon = "NO" }
                if( $stop_codon eq "YES" && $hsp->end('hit') / $hit->length > 0.85 )
                {
                    if( $orf_start < $orf_stop && abs( $orf_stop - $hsp->end('query')) < 50){ $stop_codon = "YES_REAL" }
                    if( $orf_start > $orf_stop && abs( $orf_start - $hsp->end('query')) < 50){ $stop_codon = "YES_REAL" }
                }


                my $new_desc = length($long_orf).":".$orf_start."-".$orf_stop.":".length($blast_orf).":".$hsp->start('query')."-".$hsp->end('query');
                $new_desc = $new_desc.":".$hit->name.":".$hit->length.":".$hsp->start('hit')."-".$hsp->end('hit').":".$start_codon.":".$stop_codon.":".$hsp_status;

                my $new_seq = Bio::Seq->new( -display_id => $query_id,
                                             -desc => $new_desc,
                                             -seq => $long_orf );
                $prot_file->write_seq( $new_seq );
            }
            else
            {
                my $blast_orf = $hsp->query_string;
                $blast_orf =~ s/\-//g;
                if( $blast_orf =~ /^(\w+)\*.+/ )
                {
                    if( ($hsp->end('hit') == $hit->length) && (length($blast_orf) - length($1) < 10) )
                    {
                        $blast_orf = $1;
                    }
                    else
                    {
                        next;
                    }
                }

                my $query_id = $result->query_name;
            
                my $start_codon = "NO";
                if( $hsp->start('hit') == 1 )
                {
                    if( $hsp->hit_string =~ /^M.+/ && $blast_orf =~ /^M.+/ )
                    {
                        $start_codon = "YES";
                    }
                }

                my $stop_codon = "YES";
                if( $hsp->start('query') < 4 && $result->query_length - $hsp->end('query') < 4 ){ $stop_codon = "NO" }
                if( $stop_codon eq "YES" && $hsp->end('hit') / $hit->length > 0.85 ){$stop_codon = "YES_REAL" }
                
                my ( $orf_start, $orf_stop );
                if( $hsp->strand('query') < 1 )
                {
                    $orf_start = $hsp->end('query');
                    $orf_stop = $hsp->start('query');
                }
                else
                {
                    $orf_start = $hsp->start('query');
                    $orf_stop = $hsp->end('query');
                }

                my $new_desc = length($blast_orf).":".$orf_start."-".$orf_stop.":".length($blast_orf).":".$hsp->start('query')."-".$hsp->end('query');
                $new_desc = $new_desc.":".$hit->name.":".$hit->length.":".$hsp->start('hit')."-".$hsp->end('hit').":".$start_codon.":".$stop_codon.":".$hsp_status;

                my $new_seq = Bio::Seq->new( -display_id => $query_id,
                                             -desc => $new_desc,
                                             -seq => $blast_orf );
                $prot_file->write_seq( $new_seq );

            }
        }
    }
}

