#!/usr/bin/perl -w 
#$| = 1;
use strict;
use lib ('/home/mcampbell/lib');
use Getopt::Std;
use vars qw($opt_t);
getopts('t');
use PostData;
use FileHandle;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\nscaffold_name_fixer: for fixing the ensemble scffold names
input is the maker and ensemble  GFF3 files file\n\n
this was used for a specific job but it can be tweek to do other things.             
USAGE: scaffold_name_fixer  <maker GFF3_file> <ensemble GFF3 file>\n\n";


die($usage) unless $ARGV[1];

my $A = $ARGV[0];
my $B = $ARGV[1];

my %LU;

build_LU($A);
change_name($B);

#PostData(\%LU);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub build_LU {

	my $file = shift;	
		

	my $fh = new FileHandle;
	   $fh->open($file);

	while (defined(my $line = <$fh>)){
		chomp($line);
		next unless $line =~/^scaffold/;
		last if $line =~ /^\#\#FASTA/;

		my @stuff = split(/\t/, $line);
		my ($scaff, $len) = split(/\./, $stuff[0]);
		$LU{$scaff}=$len;
	    }

	$fh->close();

}
#-----------------------------------------------------------------------------
sub change_name{

    my $file = shift;
    my $output = "Ensembl_with_exact_maker_scaffolds.gff";

    my $fh = new FileHandle;
    $fh->open($file);

    my $fho = new FileHandle;
    $fho->open(">$output");

    while (defined(my $line = <$fh>)){
	chomp($line);
	my ($scaff) = $line =~ /^(\S+)/;

	if (defined $LU{$scaff}){

		my $x = $LU{$scaff};
		my $new_line = $line;		
		$new_line =~ s/\Q$scaff/$scaff\.$x/;
		print $fho "$new_line\n";
	    }
	else {
	    print $fho "$line\n";
	}
    }
    $fh->close();
    $fho->close();
}
#-----------------------------------------------------------------------------
