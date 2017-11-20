#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_r $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u);
getopts('r:egpcmu');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
After you use fastq-dump on an sra file the reads are catted together.
This script splits the catted read in the middle (based on user supplied
read lenght -r) and prints out paired end files

sra_fastq_split_pairs.pl -r <read_length> sra_fastq_file.fastq 
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

parse($FILE);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my $out1 = $file;
    $out1 =~ s/\.fastq/_PE1\.fastq/;
    my $out2 = $file;
    $out2 =~ s/\.fastq/_PE2\.fastq/;
    my $fh = new FileHandle;
    $fh->open($file);
    my $fh1 = new FileHandle;    
    my $fh2 = new FileHandle;    
    $fh1->open(">$out1");
    $fh2->open(">$out2");
    while (defined(my $line1 = <$fh>)){
	my $line2 = <$fh>;
	my $line3 = <$fh>;
	my $line4 = <$fh>;


	chomp($line1);
	chomp($line2);
	chomp($line3);
	chomp($line4);
	$line1 =~ s/length\=\d+/length\=$opt_r/;
	$line3 =~ s/length\=\d+/length\=$opt_r/;
	
	my $r1 = substr $line2, 0, $opt_r;
	my $r2 = substr $line2, $opt_r, $opt_r;
	my $q1 = substr $line4, 0, $opt_r;
	my $q2 = substr $line4, $opt_r, $opt_r;
	#print $r1."\n";
	print $fh1 ($line1 ." 1\n");
	print $fh2 ($line1 ." 2\n");

	print $fh1 ($r1 ."\n");
	print $fh2 ($r2 ."\n");
		   
	print $fh1 ($line3 ." 1\n");
	print $fh2 ($line3 ." 2\n");
		   
	print $fh1 ($q1 ."\n");
	print $fh2 ($q2 ."\n");

    }
    $fh->close();
    $fh1->close();
    $fh2->close();
    
}
#-----------------------------------------------------------------------------

