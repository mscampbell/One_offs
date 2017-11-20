#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;
use GAL::Annotation;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
takes a  gff3 file and returns some annotation stats
\n\n";

die($usage) unless $ARGV[1];
my ( $GFF3, $FASTA ) = @ARGV;

# Make an annotation object
my $ANNOTATION = GAL::Annotation->new( $GFF3, $FASTA );

# Make a Features object using the annotation object
my $FEATURES = $ANNOTATION->features;

get_data($FEATURES);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_data{
    my $features_obj = shift;
    my $mRNA_count    =0;
    my $mRNA_length   =0;
    my $exon_count    =0;
    my $exon_length   =0;
    my $intron_count  =0;
    my $intron_length =0;
    my $tp_utr_length =0;
    my $tp_utr_count  =0;
    my $fp_utr_length =0;
    my $fp_utr_count  =0;
    my $start_count   =0;
    my $stop_count    =0;
    my $cds_length    =0;
    my $cds_count     =0;
    my $has_both      =0;
    my $single      =0;
# Do a search of the feature object and get an iterator object for all of the mRNAs
    my $mRNAs = $features_obj->search( { type => 'mRNA' } );
# Get the number of mRNAs
    $mRNA_count = $mRNAs->count;
# Go through the mRNAs one at a time
    while ( my $mRNA  = $mRNAs->next ) {
# Get the length of the mRNA
        $mRNA_length += $mRNA->length;
# Get an iterator object for the exons
        my $exons     = $mRNA->exons;
# Get the number of exons in the mRNA
        $exon_count   += $exons->count;
	$single++ if $exons->count == 1;
# Do a search of the mRNA object and get an iterator object for all of the exons
        while (my $exon   = $exons->next){
            $exon_length += $exon->length;
        }
# Even though introns do not technically exist in the gff3 file, GAL will infer them
        my $introns    = $mRNA->introns;
        $intron_count += $introns->count;
        while (my $intron   = $introns->next){
            $intron_length += $intron->length;

	}
	my $tp_utr = $mRNA->three_prime_UTR_seq; 
	if (defined($tp_utr)){
	    $tp_utr_count++; 
	    $tp_utr_length += length($tp_utr);
	}
	my $fp_utr = $mRNA->five_prime_UTR_seq;
	if (defined($fp_utr)){
	    $fp_utr_count++; 
	    $fp_utr_length += length($fp_utr);
	}
	my $cds_seq = $mRNA->CDS_seq;
	$cds_count++;
	$cds_length += length($cds_seq);
	my $start_codon = substr $cds_seq, 0, 3;
	$start_count++ if $start_codon eq 'ATG';
	my $stop_codon = substr $cds_seq, -3, 3;
	$stop_count++ if $stop_codon eq 'TAG' || $stop_codon eq 'TAA' || $stop_codon eq 'TGA';
	if ($start_codon eq 'ATG' && 
	    ($stop_codon eq 'TAG' || $stop_codon eq 'TAA' || $stop_codon eq 'TGA')){
	    $has_both++;
	}
	
#	print "$stop_codon\n";
    }
# Print out what you found
    print "\n\nfeature_type\tcount\taverage_length\tpercent_of_mRNAs_that_have_it\n";
    print "mRNA\t$mRNA_count\t". ($mRNA_length/$mRNA_count)."\n";
    print "exon\t$exon_count\t". ($exon_length/$exon_count)."\n";
    print "intron\t$intron_count\t". ($intron_length/$intron_count)."\n";
    print "3p_utr\t$tp_utr_count\t". ($tp_utr_length/$tp_utr_count)."\t". 100*($tp_utr_count/$mRNA_count)."\n";
    print "5p_utr\t$fp_utr_count\t". ($fp_utr_length/$fp_utr_count)."\t". 100*($fp_utr_count/$mRNA_count)."\n";
    print "CDS\t$cds_count\t". ($cds_length/$cds_count)."\n";
    print "start_codon\t$start_count\tNA\t". 100*($start_count/$cds_count)."\n";
    print "stop_codon\t$stop_count\tNA\t". 100*($stop_count/$cds_count)."\n";
    print "start_and_stop\t$has_both\tNA\t". 100*($has_both/$cds_count)."\n";
    print "single_exon_genes\t$single\tNA\t". 100*($single/$mRNA_count)."\n";;
}
#-----------------------------------------------------------------------------

