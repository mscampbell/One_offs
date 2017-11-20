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
input is a gff3 file and an iprscan file\n\n             
USAGE: GFF3_parser_hqgs.pl <GFF3_file> <iprscan_file>\n\n";
#The iprscan file needs to be the raw iprscan output with no header.

die($usage) unless $ARGV[1];

my $FILE    = $ARGV[0];
my $IPRFILE = $ARGV[1];

my %QUALITY;
my $GENES=0;
my $EXONS=0;
my $THQG =0;
my $SINEX=0;

parsegff3($FILE);
parseiprscan($IPRFILE);
print "NAME\t Frac_SS_confirmed\t Frac_EX_with_est_or_pro_align\t Num_pfam_dom\n";
highqset();
report();

print "\n\nThank you for choosing GFF3_parser_hqgs.pl
brought to you by the Yandell lab where 
if you\'re not happy you probably gave us 
bad data.\n\n";

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
#-----------------------------------------------------------------------------
sub parseiprscan{
    
    my $file = shift;
    
    open (my $FH, '<', $file) || die "Cant\'t find $file\n";
    
    while (defined(my $line = <$FH>)){
	chomp($line);
	loadprodom($line);
    }	
    close($FH);
}
#-----------------------------------------------------------------------------
sub loadprodom{
    my $line = shift;    
    my @ipr_col = split(/\t/, $line); 
	if ($ipr_col[9] eq "T" && ($ipr_col[8] eq "NA" || $ipr_col[8] < 1E-6)){
	    $QUALITY{genes}{$ipr_col[0]}{prodomain}++;
    }
}
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
    print "Average number of exons per gene:\t" . $EXONS/$GENES . "\n\n";
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
sub highqset{
    my @allkeys = keys %{$QUALITY{genes}};
    foreach my $name (@allkeys){
	my $x = $QUALITY{genes}{$name}{fracSSest};
	my $y = $QUALITY{genes}{$name}{fracEXestORp};
	my $z = $QUALITY{genes}{$name}{prodomain} || '';
	if ($x > 0 || $y > 0 || $z){
	    print $name;
	    if ($x < 0){
		$x = '';
	    }
	    print "\t" . $x . "\t" . $y . "\t" . $z . "\n";
	    $THQG++
	}
    }
}
#-----------------------------------------------------------------------------
