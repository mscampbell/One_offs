#!/usr/bin/perl -w 
#$| = 1;
use strict;
use lib ('/home/mcampbell/lib');
use Getopt::Std;
use vars qw($opt_t $opt_e);
getopts('te');
use PostData;
use FileHandle;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\nbinner.pl: returns values for a cumulative fraction plot based on AED or percent identity
input is a wublast output file\n
Options -t prints only the \'fraction below\' column.\n
USAGE: binner.pl  [-t] <wublast.blast or maker.gff> <bin size e.g., .025>\n
\t\tNOTE: the extentions .gff or .blast* must be present\n\n";


die($usage) unless $ARGV[1];

my $FILE = $ARGV[0];
my $BIN_UNIT = $ARGV[1];
my %HASH;
my $BIN_S = 0;
set_bins();
parse($FILE);
report() unless $opt_t;
report2() if $opt_t;
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report2{
    my $total = $HASH{'1.000000'};

    print "$FILE\n";
    print "f_below\n";
    foreach my $x (sort keys (%HASH)){
        my $frac_below = $HASH{$x}/$total;
        # more rounding to deal with floating point stuff
        my $rfrac_below = sprintf("%.3f", $frac_below);
        print "$rfrac_below\n" 
    }
}
#-----------------------------------------------------------------------------
sub report{
    my $total = $HASH{'1.000000'}; 
    print "_\t$FILE\n";
    print "f_ident\tf_below\n" if $ARGV[0] =~ /\.blast/;
    print "AED_score\tf_below\n" if $ARGV[0] =~ /\.gff/;

    foreach my $x (sort keys (%HASH)){
	my $frac_below = $HASH{$x}/$total;
        # more rounding to deal with floating point stuff
	my $rfrac_below = sprintf("%.3f", $frac_below);
	print "$x\t$rfrac_below\n";
    }
}
#-----------------------------------------------------------------------------
sub set_bins{

    while ($BIN_S <= 1){
        #I used sprintf to deal with floating point issue in perl
	my $x = sprintf("%.6f", $BIN_S);
	$HASH{$x} = 0;
	$BIN_S = $x +$BIN_UNIT;
    }
}
#-----------------------------------------------------------------------------
sub parse {
	my $file = shift;	

	my $fh = new FileHandle;
	$fh->open($file);
        # decide what kind of file it is and grab the right value
	if ($file =~ /\.blast/){
	    while (defined(my $line = <$fh>)){
		chomp($line);
                #grab the percent identity and change it
                #to a fraction before passing it to load
		my @array = split(/\t/, $line);
		load($array[10]/100);
                #grab the hsp length and change it
                #to a fraction before passing it to load
		
	    }
	}
	if ($file =~ /\.fasta/){
	    while (defined(my $line = <$fh>)){
		chomp($line);
                #grab the percent identity and change it
                #to a fraction before passing it to load
		next unless $line =~ /^>\S+\s\S+\s\S+\sAED:(\d\.\d+?)\s/;
		load($1);
                #grab the hsp length and change it
                #to a fraction before passing it to load
		
	    }
	}
        #get the AED scores out of a gff3 file	
	if ($file =~ /\.gff/){
	    while (defined(my $line = <$fh>)){
		chomp($line);
		last if $line =~ /\#\#FASTA/;
		if ($line =~ /_AED=(\d\.\d+?)\;/){
		    load($1);
		}
		else {
		    next;
		}
	    }
	}

	$fh->close();
}
#-----------------------------------------------------------------------------
sub load{
    my $value = shift;
    #go through each fo the values in the bin hash and 
    #add 1 if the value passed into the routine in less 
    #than or equal to the bin value	
    foreach my $x (keys %HASH){
	if ($value <= $x || $value eq $x){#modified after the ||
	    $HASH{$x}++;
	}
    }
}
#-----------------------------------------------------------------------------
