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
rnaseq_filter.pl scaffold_id_file gff3_file
\n\n";

my $FILE1 = $ARGV[0];
my $FILE2 = $ARGV[1];
my %LU;

die($usage) unless $ARGV[1];

parse($FILE1, 'l');
parse($FILE2, 'g');

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my $flag = shift;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	if ($flag eq 'l'){
	    $LU{$line}=1;
	}

	if ($flag eq 'g'){
	    my @stuff = split(/\t/, $line);
	    if (defined($LU{$stuff[0]}) 
		&& $stuff[1] eq 'est2genome' 
		&& $stuff[2] eq 'expressed_sequence_match'){
		my ($id) = $stuff[8] =~ /Name=(\S+)/;
		print $id ."\n";
	    }
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

