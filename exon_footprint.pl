#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\nreturns the otal numner of nucleotided covered by exons\n\n\t
exon_footprint.pl <gff3_file>
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
my $E_COUNT = 0;
my $U_COUNT = 0;
my $TOT_E = 0;
my $C_COUNT = 0;
my $TOT_C = 0;
my $TOT_U = 0;

parse($FILE);
print "\ttotal_exonic_bp\ttotal_exons\ttotal_cds_bp\ttotal_cdss\ttotal_utr_bp\ttotal_utr_exons\n";
print "$FILE\t$E_COUNT\t$TOT_E\t$C_COUNT\t$TOT_C\t$U_COUNT\t$TOT_U\n";
#print "\nexon_foot_print:$E_COUNT bp\n";
#print "utr_foot_print:$U_COUNT bp\n";
#print "total_exons:$TOT_E\n";
#print "CDS_foot_print:$C_COUNT bp\n";
#print "total_CDSs:$TOT_C\n";

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~/^\#/;
	my @array = split(/\t/, $line);

	if ($array[2] eq 'exon'){
	    my $x = $array[4] - $array[3] + 1;
	    $E_COUNT += $x;
	    $TOT_E++;
	}
	if ($array[2] eq 'five_prime_UTR' || $array[2] eq 'three_prime_UTR'){
	    my $x = $array[4] - $array[3] + 1;
	    $U_COUNT += $x;
	    $TOT_U++;
	}
	if ($array[2] eq 'CDS'){
	    my $x = $array[4] - $array[3] + 1;
	    $C_COUNT += $x;
	    $TOT_C++;
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

