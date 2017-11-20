#!/usr/bin/perl -w 
use strict;
use lib ('/sonas-hs/ware/hpc/home/mcampbel/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u $opt_a);
getopts('iegpcmua');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\t
takes a gff3 files and returns decriptive stats on the file

Total number of protein coding genes
Total number of transcripts
number of alt spliced transcripts
Average number of exons per transcript
Average number of exons per alt spliced transcript
number transcripts with 5 prime utr
Number of transcripts with 3 prime utr
Number of transcripts with AED < 0.5
\n\n";

my $FILE = $ARGV[0];
die($usage) unless $ARGV[0];
my %DATA;

parse($FILE);
report();

#PostData(\%DATA);
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub get_pc_genes{
    my %data;

    foreach my $x (keys %{$DATA{'gid'}}){
	$data{'tot_pc_genes'}++;
	$data{'tot_trans'} += $DATA{'gid'}{$x}{'num_trans'};
	if ($DATA{'gid'}{$x}{'num_trans'} > 1){
	    print $x."\n" if $opt_a;
	    $data{'tot_alt_spliced_genes'}++;
	    $data{'total_alt_spliced_trans'} += $DATA{'gid'}{$x}{'num_trans'};

	    foreach my $y (keys %{$DATA{'gid'}{$x}{'tid'}}){
		$data{'total_number_of_exons_among_alt_spliced_trans'} += 
		    $DATA{'gid'}{$x}{'tid'}{$y}{'num_exons'};
		if (defined($DATA{'gid'}{$x}{'tid'}{$y}{'num_3putr'}) && 
		    $DATA{'gid'}{$x}{'tid'}{$y}{'num_3putr'} >= 1){
		    $data{'total_number_of_alt_spliced_trans_with_3prime_utr'}++;
		}
		if (defined($DATA{'gid'}{$x}{'tid'}{$y}{'num_5putr'}) && 
		    $DATA{'gid'}{$x}{'tid'}{$y}{'num_5putr'} >= 1){
		    $data{'total_number_of_alt_spliced_trans_with_5prime_utr'}++;
		}
		$data{'total_number_of_alt_spliced_trans_with_aed_less_than_point5'}++
		    if $DATA{'gid'}{$x}{'tid'}{$y}{'aed'} < 0.5;
		$data{'total_number_of_alt_spliced_trans_with_aed_less_than_point2'}++
		    if $DATA{'gid'}{$x}{'tid'}{$y}{'aed'} < 0.2;
		$data{'total_exonic_nucleotides_alt_spliced_trans'}+=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'texon_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'texon_length'});
		$data{'total_cds_nucleotides_alt_spliced_trans'} +=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'cds_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'cds_length'});
		$data{'total_3putr_nucleotides_alt_spliced_trans'} +=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'3putr_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'3putr_length'});
		$data{'total_5putr_nucleotides_alt_spliced_trans'} +=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'5putr_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'5putr_length'});
		
		
	    }
	}
	elsif ($DATA{'gid'}{$x}{'num_trans'} == 1){
	    $data{'total_not_alt_spliced_genes'}++;
	    $data{'total_not_alt_spliced_trans'} += $DATA{'gid'}{$x}{'num_trans'};
	 
	    foreach my $y (keys %{$DATA{'gid'}{$x}{'tid'}}){
		$data{'total_number_exons_not_alt_spliced_trans'} += 
		    $DATA{'gid'}{$x}{'tid'}{$y}{'num_exons'};
		$data{'total_number_of_not_alt_spliced_trans_with_3prime_utr'}++ 
		    if (defined($DATA{'gid'}{$x}{'tid'}{$y}{'num_3putr'}) && 
			$DATA{'gid'}{$x}{'tid'}{$y}{'num_3putr'} >= 1);
		$data{'total_number_of_not_alt_spliced_trans_with_5prime_utr'}++ 
		    if (defined($DATA{'gid'}{$x}{'tid'}{$y}{'num_5putr'}) && 
			$DATA{'gid'}{$x}{'tid'}{$y}{'num_5putr'} >= 1);
		$data{'total_number_of_not_alt_spliced_trans_with_aed_less_than_point5'}++
		    if $DATA{'gid'}{$x}{'tid'}{$y}{'aed'} < 0.5;
		$data{'total_number_of_not_alt_spliced_trans_with_aed_less_than_point2'}++
		    if $DATA{'gid'}{$x}{'tid'}{$y}{'aed'} < 0.2;
		$data{'total_exonic_nucleotides_not_alt_spliced_trans'}+=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'texon_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'texon_length'});
		$data{'total_cds_nucleotides_not_alt_spliced_trans'} +=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'cds_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'cds_length'});
		$data{'total_3putr_nucleotides_not_alt_spliced_trans'} +=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'3putr_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'3putr_length'});
		$data{'total_5putr_nucleotides_not_alt_spliced_trans'} +=
		    $DATA{'gid'}{$x}{'tid'}{$y}{'5putr_length'}
		if defined($DATA{'gid'}{$x}{'tid'}{$y}{'5putr_length'});
	    }
	}
    }
#    PostData(\%data);
    return(\%data);

}
#-----------------------------------------------------------------------------
sub report{
    my ($data) = get_pc_genes();

    unless ($opt_a){

	print "Total number of protein coding genes\t$data->{'tot_pc_genes'}\n";
	print "Total number of transcripts\t$data->{'tot_trans'}\n";
	print "number of alt spliced genes\t$data->{'tot_alt_spliced_genes'}\n";
	print "number of alt spliced transcripts\t$data->{'total_alt_spliced_trans'}\n";
	print "number of not alt spliced genes\t$data->{'total_not_alt_spliced_genes'}\n";
	print "number of not alt spliced transcripts\t$data->{'total_not_alt_spliced_trans'}\n";
	print "Average number of exons per alt spliced transcript\t". $data->{'total_number_of_exons_among_alt_spliced_trans'}/$data->{'total_alt_spliced_trans'}."\n"; 
	if ($data->{'total_not_alt_spliced_trans'}){
	print "Average number of exons per not alt spliced transcript\t". $data->{'total_number_exons_not_alt_spliced_trans'}/$data->{'total_not_alt_spliced_trans'}."\n"; 
	}

	print "total number of alt spliced trans with 3 prime utr\t $data->{'total_number_of_alt_spliced_trans_with_3prime_utr'}\n";
	print "total number of not alt spliced trans with 3 prime utr\t $data->{'total_number_of_not_alt_spliced_trans_with_3prime_utr'}\n";
	print "total number of alt spliced trans with 5 prime utr\t $data->{'total_number_of_alt_spliced_trans_with_5prime_utr'}\n";
	print "total number of not alt spliced trans with 5 prime utr\t $data->{'total_number_of_not_alt_spliced_trans_with_5prime_utr'}\n";
	print "total number of alt spliced trans with AED < 0.5\t $data->{'total_number_of_alt_spliced_trans_with_aed_less_than_point5'}\n";
	print "total number of not alt spliced trans with AED < 0.5\t $data->{'total_number_of_not_alt_spliced_trans_with_aed_less_than_point5'}\n";
	print "total number of alt spliced trans with AED < 0.2\t $data->{'total_number_of_alt_spliced_trans_with_aed_less_than_point2'}\n";
	print "total number of not alt spliced trans with AED < 0.2\t $data->{'total_number_of_not_alt_spliced_trans_with_aed_less_than_point2'}\n";
	print "Average number of exonic nucleotides for alt spliced transcripts\t". $data->{'total_exonic_nucleotides_alt_spliced_trans'}/$data->{'total_alt_spliced_trans'}."\n"; 
	if ($data->{'total_not_alt_spliced_trans'}){
	print "Average number of exonic nucleotides for not alt spliced transcripts\t". $data->{'total_exonic_nucleotides_not_alt_spliced_trans'}/$data->{'total_not_alt_spliced_trans'}."\n";    
	}
	print "Average number of cds nucleotides for alt spliced transcripts\t". $data->{'total_cds_nucleotides_alt_spliced_trans'}/$data->{'total_alt_spliced_trans'}."\n"; 
	if ($data->{'total_not_alt_spliced_trans'}){
	print "Average number of cds nucleotides for not alt spliced transcripts\t". $data->{'total_cds_nucleotides_not_alt_spliced_trans'}/$data->{'total_not_alt_spliced_trans'}."\n"; 
	}
	print "Average 3 prime utr length for alt spliced transcripts\t". $data->{'total_3putr_nucleotides_alt_spliced_trans'}/$data->{'total_alt_spliced_trans'}."\n"; 
	if ($data->{'total_not_alt_spliced_trans'}){	
	print "Average 3 prime utr length for not alt spliced transcripts\t". $data->{'total_3putr_nucleotides_not_alt_spliced_trans'}/$data->{'total_not_alt_spliced_trans'}."\n"; 
	}
	print "Average 5 prime utr length for alt spliced transcripts\t". $data->{'total_5putr_nucleotides_alt_spliced_trans'}/$data->{'total_alt_spliced_trans'}."\n"; 
	if ($data->{'total_not_alt_spliced_trans'}){
	print "Average 5 prime utr length for not alt spliced transcripts\t". $data->{'total_5putr_nucleotides_not_alt_spliced_trans'}/$data->{'total_not_alt_spliced_trans'}."\n"; 
       }
    }
}
#-----------------------------------------------------------------------------
sub load_mRNAs{
    my $array_ref = shift;
    my ($gid) = $array_ref->[8] =~ /Parent=(\S+?);/;
    my ($tid) = $array_ref->[8] =~ /ID=(\S+?);/;
    my ($aed) = $array_ref->[8] =~ /_AED=(\d\.\d{2});/;
    my $mRNA_length = $array_ref->[4] - $array_ref->[3] +1;

    $DATA{'gid'}{$gid}{'num_trans'}++;
    $DATA{'gid'}{$gid}{'tid'}{$tid}{'aed'}=$aed;
    $DATA{'gid'}{$gid}{'tid'}{$tid}{'genomic_length'}=$mRNA_length;
    $DATA{'lu'}{$tid} = $gid;

}
#-----------------------------------------------------------------------------
sub load_exons{
    my $array_ref = shift;
    
    my ($tid) = $array_ref->[8] =~ /Parent=(\S+);/; #may have to adjust for ;
    my @tids = split(/,/, $tid);
    my $length = $array_ref->[4] - $array_ref->[3] +1;
    foreach my $x (@tids){
	my $gid = $DATA{'lu'}{$x};
	if ($array_ref->[2] eq 'exon'){
	    $DATA{'gid'}{$gid}{'tid'}{$x}{'num_exons'}++;
	    $DATA{'gid'}{$gid}{'tid'}{$x}{'texon_length'} += $length;
	}
	elsif ($array_ref->[2] eq 'CDS'){
	    $DATA{'gid'}{$gid}{'tid'}{$x}{'cds_length'} += $length;
	}
    }
}
#-----------------------------------------------------------------------------
sub load_UTR{
    my $array_ref = shift;
    
    my ($tid) = $array_ref->[8] =~ /Parent=(\S+);/;
    my @tids = split(/,/, $tid);
    my $utr_length = $array_ref->[4] - $array_ref->[3] +1;
    foreach my $x (@tids){
	my $gid = $DATA{'lu'}{$x};
	
	if ($array_ref->[2] eq 'five_prime_UTR'){
	    $DATA{'gid'}{$gid}{'tid'}{$x}{'num_5putr'}++;
	    $DATA{'gid'}{$gid}{'tid'}{$x}{'5putr_length'} += $utr_length;
	}
	elsif ($array_ref->[2] eq 'three_prime_UTR'){
	    $DATA{'gid'}{$gid}{'tid'}{$x}{'num_3putr'}++;
	    $DATA{'gid'}{$gid}{'tid'}{$x}{'3putr_length'} += $utr_length;
	}
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    while (defined(my $line = <$fh>)){
	chomp($line);
	my @array = split(/\t/, $line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	next unless $array[1] eq 'maker';
	load_mRNAs(\@array) if $array[2] eq 'mRNA';	
	load_exons(\@array) if $array[2] eq 'exon';
	load_UTR(\@array) if $array[2] eq 'five_prime_UTR';	
	load_UTR(\@array) if $array[2] eq 'three_prime_UTR';	
	load_exons(\@array) if $array[2] eq 'CDS';	
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

