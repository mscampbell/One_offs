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
To extract from GFF3 file source = est2genome; type = match_part. The output 
from this script will be used to extract transcripts that have 3 or more exons
for augustus training. Instead of feeding all trinity transcript to augustus, 
we only want the good ones for training, to enhance its performance. 

~/bin/MY_match-part_extrator.pl gff3_file

\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

my %DATA;

parse($FILE);
#PostData(\%DATA);
report();

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    foreach my $key (keys %DATA){
	if ($DATA{$key} > 1){
	    print "$key\n";
	}
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;
    my @gff_list;

    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);

	next if ($line !~ /^Chr/);
	last if ($line =~ /\#\#FASTA/);
	#my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $atts) = split /\t/, $line;
	@gff_list = map { split /\t/ } $line;

	if ($gff_list[2] eq 'exon'){
	    my @atts_list =  split (/;/, $gff_list[8]);
	    foreach my $i (@atts_list) {
		if ( $i =~ /Parent=(\S+)/) {
		    $DATA{$1}++;
		}	
	    }		
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

