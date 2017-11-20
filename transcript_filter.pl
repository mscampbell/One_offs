#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');
use PostData;
use FileHandle;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\ttranscript_filter.pl: filters the transcripts in a gff3 file 
based on a list of ids.\n
USAGE: star_keeper.pl <id_file> <gff3_file>\n";

die($usage) unless $ARGV[1];

my $IDS = $ARGV[0];
my $GFF = $ARGV[1];
my %LU;
my %LU2;
load_lu($IDS); 
load_lu2($GFF); 

parse($GFF);

#PostData(\%LU);
#PostData(\%LU2);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub load_lu2{
    my $file = shift;       

    my $fh = new FileHandle;
    $fh->open($file);

    while (defined(my $line = <$fh>)){
	chomp($line);

	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	my @array = split(/\t/, $line);

	if ($array[2] =~ 'mRNA'){
	    my ($tid) = $line =~ /ID=(.+?);/;
	    my ($gid) = $line =~ /Parent=(.+?);/;
	    if ($LU{'tid'}{$tid}){
		$LU2{'tid'}{$tid} = $gid;
		$LU2{'gid'}{$gid} = $tid;
	    }
	    else {next;}
	}
    }
}
#-----------------------------------------------------------------------------
sub load_lu{
    my $file = shift;       

    my $fh = new FileHandle;
    $fh->open($file);

    while (defined(my $line = <$fh>)){
	chomp($line);
	my ($gid) = $line =~ /(.+?)\.\d/;
	$LU{'tid'}{$line} = $gid;
	$LU{'gid'}{$gid} = $line;

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
	
	if ($array[2] =~ 'gene'){
	 my @kvp =split(/\;/, $array[8]);
            my @id = split(/\=/, $kvp[0]);
	 if ($LU{'gid'}{$id[1]} && $LU2{'gid'}{$id[1]}){
            print $line."\n";
	    }
            else {next};    
	}
    	elsif ($array[2] =~ 'mRNA'){
	    my @kvp =split(/\;/, $array[8]);
	    my @id = split(/\=/, $kvp[0]);
	    if ($LU{'tid'}{$id[1]} && $LU2{'tid'}{$id[1]}){
	    print $line."\n";
	    }
	    else {next};
	}
	elsif ($array[2] =~ 'exon' || $array[2] =~ 'CDS'){
	    my @kvp =split(/\;/, $array[8]);
            my @id = split(/\=/, $kvp[1]);
	    my @list = split(/\,/, $id[1]);
	    foreach my $x (@list){
		if ($LU{'tid'}{$x} && $LU2{'tid'}{$x}){
		    print $line."\n";
		}
		else {next};
	    }
	}
	elsif ($array[2] =~ 'five_prime_UTR' || $array[2] =~ 'three_prime_UTR'){
	    my @kvp =split(/\;/, $array[8]);
            my @id = split(/\=/, $kvp[1]);
	    my @list = split(/\,/, $id[1]);
	    foreach my $x (@list){
		if ($LU{'tid'}{$x} && $LU2{'tid'}{$x}){
		    print $line."\n";
		}
		else {next};
	    }
	}
	else {
	   print $line."\n";
	}	
    }
}
#-----------------------------------------------------------------------------
