#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');
use PostData;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\nfinds the scaffols witht the most genes\n\n
USAGE: dense_scaffold.pl <gff file>\n\n";             


die($usage) unless $ARGV[0];

my $FILE1 = $ARGV[0];
my %DATA;

parse($FILE1);
report();


#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub parse {
    
    my $file = shift;
    
    open (my $FH, '<', $file) || die "Cant\'t find $file\n";
    
    while (defined(my $line = <$FH>)){
        chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
        next if $line =~ /^\s*$/;    
	
	my @stuff = split(/\t/,$line);
#	# make the maker scaffold ids look like the ones in ensembl
	my @idline = split(/\./, $stuff[0]);	
foreach my $e ($stuff[2]){
	    if ($e =~ "mRNA"){
		$DATA{$idline[0]}++
		}
	    
	}
    }
    
    close($FH);
}

#-----------------------------------------------------------------------------
sub report{
    foreach my $s (sort { $DATA{$b} <=> $DATA{$a} }keys (%DATA)){
	my $x = $DATA{$s};
	# the number below can be adjusted to catch scaffolds with fewer genes
	if ($x > 30){
	    print $s ."\t". $x . "\n";
	}
    }   
}


#-----------------------------------------------------------------------------
    
