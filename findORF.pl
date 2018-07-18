#!/usr/bin/perl -w

use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::SeqIO;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

my $file = $ARGV[0];
my $report = Bio::SearchIO->new( -format => 'blast', -file => $ARGV[0] );
my $db = Bio::DB::Fasta->new( $ARGV[1] );
my ($counter, $hit_counter, $one_hsp_counter, $incomplete_counter, $probably_complete, $def_complete);
my ($unique_def_complete, $unique_probably_complete, $retro_counter, $potential_complete);
my $hairpin_counter;
my %prots = ();
my %genes = ();
my $gene_counter;
my %mito_genes = ();
my ($mito_counter_seqs, $mito_counter_genes);
my %chloro_genes = ();
my ($chloro_counter_seqs, $chloro_counter_genes);
my %orf_genes = ();
my $orf_gene_counter;
open( OUTTY, ">".$ARGV[0].".orf" );

while( my $result = $report->next_result )
{
    $counter++;
    if( $result->num_hits > 0 )
    {
        $hit_counter++;
        my $hit = $result->next_hit;
        my $gene_name = "";
        if( $hit->description =~ /.+gene:(\S+)\s+.*/ )
        {
            $gene_name = $1;
            if( !defined( $genes{$gene_name} ))
            {
                $genes{$gene_name} = 1;
                $gene_counter++;
            }
        }
        
        if( $hit->description =~ /.+chromosome:\S+:Mt:\d+.+/ ||
            $hit->description =~ /.+chromosome:\S+:mitochondrion:\d+.+/ )
        {
            $mito_counter_seqs++;
            if( !defined( $mito_genes{$gene_name} ))
            {
                $mito_genes{$gene_name} = 1;
                $mito_counter_genes++;
            }
        }

        if( $hit->description =~ /.+chromosome:\S+:Pt:\d+.+/ ||
            $hit->description =~ /.+chromosome:\S+:chloroplast:\d+.+/ )
        {
            $chloro_counter_seqs++;
            if( !defined( $chloro_genes{$gene_name} ))
            {
                $chloro_genes{$gene_name} = 1;
                $chloro_counter_genes++;
            }
        }

        if( $hit->description =~ /transposon/ || 
            $hit->description =~ /transposable/ ){ $retro_counter++ }

        if( $hit->num_hsps == 1 )
        {
            my $hsp = $hit->next_hsp;
            $one_hsp_counter++;
            if( $hsp->start('query') < 4 || 
                $hsp->start('query') > ($result->query_length - 4 ) ||
                $hsp->end('query') < 4 ||
                $hsp->end('query') > ($result->query_length - 4 ) )
            {
                $incomplete_counter++;
            }
            else
            {
                $potential_complete++;
                my $seq = $db->get_Seq_by_id( $result->query_name );
                my $temp_name = $file."_ORF_".$potential_complete;
                my $temp_file = Bio::SeqIO->new( -file => ">".$temp_name, -format => 'fasta' );
                $temp_file->write_seq( $seq );
                my $frame = $hsp->query->frame;
                my $strand = $hsp->strand;
                $frame++;
                my $frame_name = "frame ".$frame;
                if( $strand < 1 )
                {
                    if( $frame == 1 ){ $frame_name = "frame 5" }
                    if( $frame == 2 ){ $frame_name = "frame 4" }
                    if( $frame == 3 ){ $frame_name = "frame 6" }
                }
                my $cmd = "sixpack -mstart -auto -sequence $temp_name -outfile $temp_name.out -outseq $temp_name.orf";
                system( $cmd );
                my @seqs = `grep "$frame_name" $temp_name.orf`;
                my $orf_name = "none";
                my $orf_length = 1;
                my $unique = 0; #false
                foreach my $orf_line (@seqs)
                {
                    chomp $orf_line;
                    $orf_line =~ /^\>(\S+)\s+.+threshold\s\d+\,\s(\d+)aa.*/;
                    if( $2 > $orf_length )
                    {
                        $orf_name = $1; 
                        $orf_length = $2;
                    }
                }
                if( $orf_name ne "none" )
                {
                    my $orf_db = Bio::DB::Fasta->new( $temp_name.".orf" );
                    my $orf_seq = $orf_db->get_Seq_by_id( $orf_name );
                    my $orf_prot = $orf_seq->seq;
                    if( !defined( $prots{$orf_prot} ))
                    {
                        $unique = 1;
                        $prots{$orf_prot} = 1;
                    }
                }
                $cmd = "rm $temp_name*";
                system( $cmd );
                my $hit_length = $hit->length;
                if( (max $orf_length, $hit_length) - (min $orf_length, $hit_length) < 20 )
                {
                    $def_complete++;
                    if( $unique == 1 )
                    { 
                        $unique_def_complete++;
                        if( !defined( $orf_genes{$gene_name} ))
                        {
                            $orf_genes{$gene_name} = 1;
                            $orf_gene_counter++;
                        }
                    }
                }
                elsif( (min $orf_length, $hit_length) / (max $orf_length, $hit_length) * 100 > 75 )
                {
                    $probably_complete++;
                    if( $unique == 1 )
                    { 
                        $unique_probably_complete++;
                        if( !defined( $orf_genes{$gene_name} ))
                        {
                            $orf_genes{$gene_name} = 1;
                            $orf_gene_counter++;
                        }
                    }
                }
            }
        }
        else
        {
            my $frame_checker = 0;
            while( my $hsp = $hit->next_hsp )
            {
                $frame_checker += $hsp->strand;
            }
            if( $frame_checker != 0 && $frame_checker != $hit->num_hsps ){ $hairpin_counter++ }
        }
    }    
}

print OUTTY "\n\nCounter: $counter\nHit Counter: $hit_counter\nMatching Genes: $gene_counter\n";
print OUTTY "$mito_counter_seqs transcripts matching to $mito_counter_genes mitochondrial genes\n";
print OUTTY "$chloro_counter_seqs transcripts matching to $chloro_counter_genes chloroplast genes\n";
print OUTTY "Transposons/Retrotransposons: $retro_counter\n";
print OUTTY "Potential hairpins/palindromes/misassembly: $hairpin_counter\n";
print OUTTY "One HSP Counter: $one_hsp_counter\n";
print OUTTY "Incomplete Counter: $incomplete_counter\n";
print OUTTY "Probably Complete: $probably_complete\nUnique Probably Complete: $unique_probably_complete\n";
print OUTTY "Definately Complete: $def_complete\nUnique Definately Complete: $unique_def_complete\n";
print OUTTY "Unique Potential/Definates matching to $orf_gene_counter genes\n";
