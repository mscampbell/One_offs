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
my $usage = "\n\n\t
I used this to get one to one overlapping genes in arabidopsis
\n\n";
my %CB;
my %CA;
my %B2A;
my %A2B;

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

parse($FILE);
something();
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub something{
    foreach my $x (keys %CB){
	print $B2A{$x}."\n" if $CB{$x} == 1;
	print $x."\n" if $CB{$x} == 1;
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	my @stuff = split(/\t/, $line);
	$CB{$stuff[1]}++;
	$B2A{$stuff[1]}=$stuff[0];
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

