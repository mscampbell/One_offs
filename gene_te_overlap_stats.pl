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
Parses bedtool intercest output file from two gff3 files

example bedtools command: 
\tbedtools intersect -wo -a B73.v4.DEF.cleaned.rep.gff -b B73.Mhap2.quiver.tirvish.tir11.gff3

\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];

my $data = parse($FILE);
report($data);
#PostData($data);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub report{
    print "#gene_id\tnumber_of_exons_overlaped\tCDS_overlap\tids_of_te_features_overlapping_the_gene\tfraction_of_the_gene_overlapped\n";
    my $data = shift;
    foreach my $gene (keys %$data){
	my $number_of_exons_overlaped = 0;
	my $cds_overlap = 0;
	my @ids_of_te_features_overlapping_the_gene;
	my @fraction_of_the_gene_overlapped;
	foreach my $type (keys %{$data->{$gene}}){
	    if ($type eq 'exon'){
		foreach my $ex_id (keys %{$data->{$gene}->{$type}}){
		    $number_of_exons_overlaped++
		}
	    }
	    elsif ($type eq 'CDS'){
		$cds_overlap = 1;
	    }
	    elsif ($type eq 'gene'){
		foreach my $te_id (keys %{$data->{$gene}->{$type}->{$gene}->{te_id}}){
		    my $gene_len = $data->{$gene}->{$type}->{$gene}->{'f_len'};
		    my $te_overlap = $data->{$gene}->{$type}->{$gene}->{te_id}->{$te_id}->{$te_id};
		    my $frac_overlap = $te_overlap/$gene_len;
		    push (@ids_of_te_features_overlapping_the_gene, $te_id);
		    push (@fraction_of_the_gene_overlapped, $frac_overlap);
		}
	    }
	}
	print "$gene\t$number_of_exons_overlaped\t$cds_overlap\t";
	print join(",", @ids_of_te_features_overlapping_the_gene);
	print "\t";
	print join(",", @fraction_of_the_gene_overlapped);

	print "\n";

	
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my %data;
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	my @col = split("\t", $line);
	my $g_type = $col[2];
	my ($gf_id) = $col[8] =~ /ID=(\S+?);/;
	my $gid;
	if ($gf_id =~ /-mRNA.*/){
	#    ($gf_id) = $gid =~ /(\S+?)/;

	$gid = $gf_id;
	$gid =~ s/-mRNA.*//;
	}
	else{$gid = $gf_id;}
	
	my $te_type = $col[11];
	my $te_id;
	my $tef_id;
	next if $col[17] !~ /ID=/; 
	if ($col[17] =~ /Parent=(\S+?);/){
	    $te_id = $1;
	    ($tef_id)= $col[17] =~ /ID=(\S+?);/;  
	} 
	elsif($col[17] =~ /ID=(\S+)/){
	    $te_id = $1;
	    $tef_id = $te_id;
	}
	else{print "saw this $col[17]\n";}

	my $len_ovl = $col[18];
	my $len_gf  = $col[4] - $col[3] + 1;

#	print "type\t$g_type\ngid\t$gid\nte_type\t$te_type\nte_id\t$te_id\n\n";
	
	$data{$gid}{$g_type}{$gf_id}{'f_len'} = $len_gf;
	$data{$gid}{$g_type}{$gf_id}{'te_id'}{$te_id}{$tef_id} += $len_ovl;

    }
    $fh->close();
    return \%data;
}
#-----------------------------------------------------------------------------

