#!/usr/bin/perl -w 
use strict;
use Getopt::Std;
use vars qw($opt_e);
getopts('e');

use warnings;
use lib ('/Users/mcampbell/mcampbell/lib');
use PostData;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "\n\nGFF_parser: outputs a HQ gene set 
input is a gff3 file\n\n             
USAGE: GFF3_parser_hqgs.pl <GFF3_file>\n\n";


die($usage) unless $ARGV[0];

my $FILE    = $ARGV[0];


my %QUALITY;
my $GENES=0;
my $EXONS=0;
my $THQG =0;
my $SINEX=0;
my $PROD =0;

parsegff3($FILE);
#print "NAME\tFrac_SS_confirmed\tFrac_EX_with_est_or_pro_align\n";
#highqset();
report();

print "\n\nThank you for choosing GFF3_parser_hqgs.pl
brought to you by the Yandell lab where 
if you\'re not happy you probably gave us 
bad data.\n\n";
#PostData(\%QUALITY);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parsegff3 {
    
    my $file = shift;
    
    open (my $FH, '<', $file) || die "Cant\'t find $file\n";
    
    while (defined(my $line = <$FH>)){
	chomp($line);
	next if $line =~/^#|^>|^[A-Z]/;
	
	genecount($line);
	loadquality($line);
    }	
    close($FH);
}
#gff3 files are 9 columns
#c1=seqid (landmark ID, usualy the scaffold or contig name in MAKER) 
#c2=source (the thing that generated the feature)
#c3=type (SO number)
#c4=start (specific to the scaffold or contig)
#c5=end (specific to the scaffold or contig)
#c6=score (e.g., this is and E value for blast)
#c7=strand (+,-,or ?)
#c8=phase (reading frame 1,2,or 3)
#c9=attributes (gene names and maker generated quality scores are here)
#-----------------------------------------------------------------------------
sub genecount{
    my $line = shift;    
    my @GFF_col = split(/\t/, $line); 
    $GENES++ if $GFF_col[2] =~/gene/;
    $EXONS++ if $GFF_col[2] =~/exon/;
}
#-----------------------------------------------------------------------------
sub report{
    print "\n\nTotal high quality genes:\t" . $THQG . "\n";
    print "Total genes:\t" . $GENES . "\n";     
    print "Total exons:\t" . $EXONS . "\n";     
    print "Average number of exons per gene:\t" . $EXONS/$GENES . "\n";
}
#-----------------------------------------------------------------------------
sub loadquality{
    my $line = shift;    
    my @GFF_col = split(/\t/, $line); 
    
    if ($GFF_col[2] eq 'mRNA'){ 
	my @mrna9 = split(/\;/, $GFF_col[8]); #@mrana has 6 elements 
	$mrna9[0] =~ s/ID=//;
	$mrna9[1] =~ s/Parent=//;
	$mrna9[2] =~ s/Name=//;
	$mrna9[3] =~ s/_AED=//;
	$mrna9[4] =~ s/_eAED=//;
	$mrna9[5] =~ s/_QI=//;
	my @qi = split(/\|/, $mrna9[5]);
	my @name = split(/-/, $mrna9[2]);
	my $scaffold = $name[1];
	$QUALITY{genes}{$mrna9[0]}{scaffold}=$scaffold;
	$QUALITY{genes}{$mrna9[0]}{lenFPutr}=$qi[0];
	$QUALITY{genes}{$mrna9[0]}{fracSSest}=$qi[1];
	$QUALITY{genes}{$mrna9[0]}{fracEXest}=$qi[2];
	$QUALITY{genes}{$mrna9[0]}{fracEXestORp}=$qi[3];
	$QUALITY{genes}{$mrna9[0]}{fracSSgp}=$qi[4];
	$QUALITY{genes}{$mrna9[0]}{fracEXgp}=$qi[5];
	$QUALITY{genes}{$mrna9[0]}{numEX}=$qi[6];
	$QUALITY{genes}{$mrna9[0]}{lenTPutr}=$qi[7];
	$QUALITY{genes}{$mrna9[0]}{lenPRO}=$qi[8];
    }
    
}
#-----------------------------------------------------------------------------
##########left out because the lamprey GFF has a different column 9 order#####
#sub highqset{
#    my @allkeys = keys %{$QUALITY{genes}};
#    foreach my $name (@allkeys){
#	my $x = $QUALITY{genes}{$name}{fracSSest};
#	my $y = $QUALITY{genes}{$name}{fracEXestORp};	
#	if ($x > 0 || $y > 0){
#	    print $name;
#	    if ($x < 0){
#		$x = '';
#	    }
#	    print "\t" . $x . "\t" . $y;
#	    print "\n";
#	    $THQG++
#	}
#    }
#}
#-----------------------------------------------------------------------------
	
