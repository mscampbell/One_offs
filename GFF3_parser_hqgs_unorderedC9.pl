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
USAGE: GFF3_parser_hqgs_gffOnly.pl <GFF3_file>\n\n";
#The iprscan file needs to be the raw iprscan output with no header.

die($usage) unless $ARGV[0];

my $FILE    = $ARGV[0];


my %QUALITY;
my $GENES=0;
my $EXONS=0;
my $MRNA =0;
my $THQG =0;
my $SINEX=0;
my $PROD =0;

parsegff3($FILE);
print "NAME\tFrac_SS_confirmed\tFrac_EX_with_est_or_pro_align\n";
highqset();
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
    $MRNA++ if $GFF_col[2] =~/mRNA/;
}
#-----------------------------------------------------------------------------
sub report{
    print "\n\nTotal high quality genes:\t" . $THQG . "\n";
    print "Total genes:\t" . $GENES . "\n";     
    print "Total transcripts:\t" . $MRNA . "\n";
    print "Total exons:\t" . $EXONS . "\n";     
    print "Average number of exons per transcript:\t" . $EXONS/$MRNA . "\n";

}
#-----------------------------------------------------------------------------
sub loadquality{
    my $line = shift;    
    my @GFF_col = split(/\t/, $line); 
    my $gname;
    my $scaffold;
    
    if ($GFF_col[2] eq 'mRNA'){ 
	my @mrna9 = split(/\;/, $GFF_col[8]); #@mrana9 has 6 elements 
	
	foreach my $att (@mrna9){
	    if ($att =~ /^ID=/){
		$att  =~ s/ID=//;
		$gname = $att;
		my @lname = split(/-/, $att);
		$scaffold = $lname[1];    
	    }
	    elsif ($att =~ /_QI=/){
		$att =~ s/_QI=//;
		my @qi = split(/\|/, $att);
		
		$QUALITY{genes}{$gname}{scaffold}=$scaffold;
		$QUALITY{genes}{$gname}{lenFPutr}=$qi[0];
		$QUALITY{genes}{$gname}{fracSSest}=$qi[1];
		$QUALITY{genes}{$gname}{fracEXest}=$qi[2];
		$QUALITY{genes}{$gname}{fracEXestORp}=$qi[3];
		$QUALITY{genes}{$gname}{fracSSgp}=$qi[4];
		$QUALITY{genes}{$gname}{fracEXgp}=$qi[5];
		$QUALITY{genes}{$gname}{numEX}=$qi[6];
		$QUALITY{genes}{$gname}{lenTPutr}=$qi[7];
		$QUALITY{genes}{$gname}{lenPRO}=$qi[8];
	    }
	    elsif ($att =~ /_AED=/){
		$att =~ s/_AED=//;
		$QUALITY{genes}{$gname}{AED}=$att;
	    }

	    elsif ($att =~ /Ontology_term=GO/){
		$att =~ s/Ontology_term=GO//;
		$QUALITY{genes}{$gname}{GoTerms}=$att;
	    }
	}
    }
}
#-----------------------------------------------------------------------------
sub highqset{
    my %hqt;    
    my @allkeys = keys %{$QUALITY{genes}};
    foreach my $name (@allkeys){
	my $x = $QUALITY{genes}{$name}{fracSSest};
	my $y = $QUALITY{genes}{$name}{fracEXestORp};
	
	if ($x > 0 || $y > 0){
	    print $name;
	    if ($x < 0){
		$x = '';
	    }
	    print "\t" . $x . "\t" . $y;
	    print "\n";
	    $THQG++;
	    $hqt{$name}++;
	}
    }
    my @keys = keys(%hqt);
    foreach my $key (@keys){
	if ($hqt{$key} > 1){
	    print $key
	}
	
	else{ 
	    print '';
	}
    }
}
#-----------------------------------------------------------------------------

