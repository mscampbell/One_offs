#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('iegpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
It looks like this script takes a gff and changes the ID of each feature to 
reflect the file it came from. I probably used this for somthing I needed
to see in a browser.
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

parse($FILE);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	last if $line =~ /\#\#FASTA/;
	next if $line =~ /\#/;
	
	my @col = split(/\t/, $line);
	if ($col[1] eq 'maker'){
	    $line =~ s/ID=/ID=$file/ if $line =~ /ID=/;
	    $line =~ s/Parent=/Parent=$file/ if $line =~ /Parent=/;
	    $line =~ s/Name=/Name=$file/ if $line =~ /Name=/;
	    print $line."\n";
	}
	else {print $line."\n" unless $opt_g;}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

