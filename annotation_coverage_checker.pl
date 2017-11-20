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
my $usage = "\n\nannotation_coverage_checker.pl: calculates how much of a 
subject sequence is overlaped by the best hit query sequence.  
input is a tab delimited wublast report\n\n
USAGE: annotation_coverage_checker.pl  <.blast file>\n\n";


die($usage) unless $ARGV[0];

my $A = $ARGV[0];

my %DATA;

build_DATA($A);

#PostData(\%DATA);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub build_DATA {

	my $file = shift;	
		

	my $fh = new FileHandle;
	   $fh->open($file);

	while (defined(my $line = <$fh>)){
		chomp($line);
		
		my @stuff = split(/\t/, $line);
		my $alignlen = $stuff[6];
		my $ggaplen = $stuff[13];
		my $fracc = ($alignlen - $ggaplen)/$alignlen;
		my $covered = ($alignlen - $ggaplen);
		$DATA{$stuff[0]}{covered} = $covered;
		$DATA{total_bases_covered} += $covered;
		$DATA{total_bases} += $alignlen;
		push (@{$DATA{fracc_array}}, $fracc);
		print "$stuff[17]\n";
	    }

	$fh->close();

}
##-----------------------------------------------------------------------------
#sub change_name{
#
#    my $file = shift;
#    my $output = "Ensembl_with_exact_maker_scaffolds.gff";
#
#    my $fh = new FileHandle;
#    $fh->open($file);
#
#    my $fho = new FileHandle;
#    $fho->open(">$output");
#
#    while (defined(my $line = <$fh>)){
#	chomp($line);
#	my ($scaff) = $line =~ /^(\S+)/;
#
#	if (defined $LU{$scaff}){
#
#		my $x = $LU{$scaff};
#		my $new_line = $line;		
#		$new_line =~ s/\Q$scaff/$scaff\.$x/;
#		print $fho "$new_line\n";
#	    }
#	else {
#	    print $fho "$line\n";
#	}
#    }
#    $fh->close();
#    $fho->close();
#}
##-----------------------------------------------------------------------------
