#!/usr/bin/perl  

# Use the module in GAL that lets you make the annotation object
use GAL::Annotation;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
This is to get you started using GAL. This script will return
the total number and average lengths of specified features in
a gff3 file
Simple_GAL_script.pl GFF3.gff FASTA.fa
\n\n";

die($usage) unless $ARGV[1];
my ( $GFF3, $FASTA ) = @ARGV;

# Make an annotation object
my $ANNOTATION = GAL::Annotation->new( $GFF3, $FASTA );

# Make a Features object using the annotation object
my $FEATURES = $ANNOTATION->features;

get_length_and_count($FEATURES);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_length_and_count{
    my $features_obj = shift;
    my $mRNA_count    =0;
    my $mRNA_length   =0;
    my $exon_count    =0;
    my $exon_length   =0;
    my $intron_count  =0;
    my $intron_length =0;

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
    }

# Print out what you found
    print "\n\nfeature type\tcount\taverage_length\n";
    print "mRNA\t$mRNA_count\t". int($mRNA_length/$mRNA_count)."\n";
    print "exon\t$exon_count\t". int($exon_length/$exon_count)."\n";
    print "intron\t$intron_count\t". int($intron_length/$intron_count)."\n";
}
