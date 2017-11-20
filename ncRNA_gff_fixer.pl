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
my $usage = "\n\ntakes Yanni\'s ncRNA gff3 file and makes it look like 
protein coding gene annotations.
ncRNA_gff_fixer.pl <gff3_file>

\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

parse($FILE);

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub make_gene{
    my $array = shift;
#    print $array->[0]."\n";
    my $c1 = $array->[0];
    my $c2 = $array->[1];
    my $c3 = 'gene';
    my $c4 = $array->[3];
    my $c5 = $array->[4];
    my $c6 = $array->[5];
    my $c7 = $array->[6];
    my $c8 = $array->[7];
    my $c9 = 'ID='."$array->[2]-$array->[0]-$array->[3]-$array->[4];$array->[8]";
    my @line = ($c1,$c2,$c3,$c4,$c5,$c6,$c7,$c8,$c9);
    return \@line;
}
#-----------------------------------------------------------------------------
sub make_mrna{
    my $array = shift;
#    print $array->[0]."\n";                                                                             
    my $c1 = $array->[0];
    my $c2 = $array->[1];
    my $c3 = 'mRNA';
    my $c4 = $array->[3];
    my $c5 = $array->[4];
    my $c6 = $array->[5];
    my $c7 = $array->[6];
    my $c8 = $array->[7];
    my ($parent) = $array->[8] =~ /ID=(.+?);/;
    my ($x) = $parent =~ /(.+?)-/;
    my $id = 'ID='."$x-$array->[0]-$array->[3]-$array->[4]-mRNA";
    my $c9 = "$id;Parent=$parent";
    my @line = ($c1,$c2,$c3,$c4,$c5,$c6,$c7,$c8,$c9);
    return \@line;
}
#-----------------------------------------------------------------------------
sub make_exon{
    my $array = shift;
#    print $array->[0]."\n";                                                                           
                                                                                                        
    my $c1 = $array->[0];
    my $c2 = $array->[1];
    my $c3 = 'exon';
    my $c4 = $array->[3];
    my $c5 = $array->[4];
    my $c6 = $array->[5];
    my $c7 = $array->[6];
    my $c8 = $array->[7];
    my ($parent) = $array->[8] =~ /ID=(.+?);/;
    my ($x) = $parent =~ /(.+?)-/;
    my $id = 'ID='."$x-$array->[0]-$array->[3]-$array->[4]-exon";
    my $c9 = "$id;Parent=$parent";
    my @line = ($c1,$c2,$c3,$c4,$c5,$c6,$c7,$c8,$c9);
    return \@line;
}

#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	next if $line =~ /^\#/;
	last if $line =~ /^\#\#FASTA/;
	my @array = split(/\t/, $line);
	my $gene_line = make_gene(\@array);
	my $mrna_line = make_mrna($gene_line);
	my $exon_line = make_exon($mrna_line);
	print join("\t", @$gene_line) ."\n";
	print join("\t", @$mrna_line) ."\n";
	print join("\t", @$exon_line) ."\n";

    }
    $fh->close();
}
#-----------------------------------------------------------------------------

