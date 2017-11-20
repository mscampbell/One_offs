#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_n);
getopts('n:');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\test_filter_for_gene_finder_training.pl: takes a
\tMAKER generated gff3 file and returns the ids of 
\tthe ests with a given number of \'exons\'
\n\nest_filter_for_gene_finder_training.pl -n <min number of exons> MAKER.gff
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
die($usage) unless $opt_n;

my %DATA;

parse($FILE);
#PostData(\%DATA);
report();
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{

    foreach my $x (keys %DATA){
	foreach my $y (keys %{$DATA{$x}}){
	print "$x\t$DATA{$x}{$y}\n" if $DATA{$x}{$y} >= $opt_n;
	
	}
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @array = split(/\t/, $line);
	next unless $array[1] eq 'est2genome' 
	            && $array[2] eq 'match_part';

	my ($parent) = $array[8] =~ /Parent=(\S+?);/; 
	my ($name) = $array[8] =~ /Target=(\S+)/; 

	$DATA{$name}{$parent}++;

    }
    $fh->close();
}
#-----------------------------------------------------------------------------

