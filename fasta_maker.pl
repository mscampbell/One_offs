#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');
use PostData;
use Getopt::Std;
use vars qw($opt_i $opt_e $opt_g $opt_p $opt_c $opt_m $opt_u $opt_s);
getopts('iegpcmus:');
use FileHandle;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "\n\n\tfasta_maker.pl: generates fasta sequences given a maker
\tgenerated gff3 file and a multi fasta file of the assembly\n\n
USAGE: fasta_maker.pl [igpemu] <gff3_file> <multi_fasta_file>

OPTIONS: -i returns intronic sequences in fasta format
         -g returns the whole gene in fasta format
         -p returns intergenic sequences in fasta format
         -e rerutns exonic sequences in fasta format
         -m returns the pre_spliced mRNA in fasta format
         -u returns the UTRs in fasta format
         -s<string> the first two or three characters of the scaffold names 

EXAMPLE: fasta_maker.pl -p -s Chr test_genes.gff TAIR10_chr2.fa > output.fa
\n\n";
die($usage) unless $ARGV[1];
die "\nCDS is broken\n" if $opt_c;

my $GFF   = $ARGV[0]; 
my $FASTA = $ARGV[1];
my %SEQ;
my $DEF;
my %ST_ST;
my $total_lenght = 0;
my %P_GENES;
parse($GFF, 'g');
parse($FASTA, 'f');

print_fasta('g') if $opt_g;
print_fasta('e') if $opt_e;
print_fasta('m') if $opt_m;
print_fasta('i') if $opt_i;
print_fasta('u') if $opt_u;
print_fasta('c') if $opt_c;
print_ig('p')    if $opt_p;



 
#PostData(\%SEQ);
print STDERR "\ntotal_nucleotides:$total_lenght\n";
#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
sub reverse_comp{

    my $seq = shift;
    my $r_seq;

    $r_seq = reverse($seq);
    $r_seq =~ tr/ACGTYRKMB/TGCARUMKV/;

    return $r_seq;
}
#-----------------------------------------------------------------------------
sub get_feature{

    my $flag = shift;
    my $feature;

    if ($flag eq 'g'){
        $feature = 'whole_gene';
    }
    if ($flag eq 'm'){
        $feature = 'whole_mRNA';
    }
    if ($flag eq 'e'){
        $feature = 'all_exons_and_utrs';
    }
    if ($flag eq 'five_prime_UTR'){
        $feature = 'all_utrs';
    }
    if ($flag eq 'three_prime_UTR'){
        $feature = 'all_utrs';
    }
    if ($flag eq 'i'){
        $feature = 'all_introns';
    }
    if ($flag eq 'c'){
        $feature = 'all_CDSs';
    }
    if ($flag eq 'p'){
        $feature = 'intergenic';
    }
    return $feature;
}
#-----------------------------------------------------------------------------
sub print_it{

    my $id    = shift;
    my $start  = shift;
    my $stop   = shift;
    my $scaf   = shift;
    my $strand = shift;
    my $flag   = shift;
    my $count = shift;
    my $feature=get_feature($flag);
    
    if ($opt_i || $opt_p){
	$start++;
	$stop--;
    }

    my $length = $stop - $start +1;
    $total_lenght += $length;
    my $index  = $start -1;

    if ($strand eq '-'){
	print ">$id count_$count B:$start E:$stop $feature\t-_strand_feature_rev_comp\n";
	my $seq = substr($SEQ{$scaf}, $index, $length);
	my $rev_comp = reverse_comp($seq);
	print $rev_comp."\n"
    }
    if ($strand eq '+'){
	print ">$id count_$count B:$start E:$stop $feature\t\+_strand_feature\n";
	print substr($SEQ{$scaf}, $index, $length)."\n"; 
    }
}
#-----------------------------------------------------------------------------
sub print_ig{

    my $flag  = shift;
    my $count = 0;
    my @array;
    my $id;
    foreach my $scaf (keys %P_GENES){
	@array = sort {$a <=> $b} @{$P_GENES{$scaf}};
	shift @array;
	while ($array[1]){
	    my $start = shift(@array);
	    my $stop  = shift(@array);
	    $count++;
	    $id = $scaf."_".$start."_".$stop;
	    print_it($id, $start, $stop, $scaf, '+', $flag, $count);
	}
    }
}
#-----------------------------------------------------------------------------
sub print_fasta{

    my $flag = shift;
    my $count = 0;
    my @array;
    
    foreach my $scaf (keys %ST_ST){
	foreach my $id (keys %{$ST_ST{$scaf}{$flag}}){
	    @array = sort {$a <=> $b} @{$ST_ST{$scaf}{$flag}{$id}{'stst'}};
	    if ($flag eq 'i'){
		shift @array;
		}
	    
	    while ($array[1]){
		my $start = shift(@array);
		my $stop  = shift(@array);
		my $strand = $ST_ST{$scaf}{'strand'}{$id};
		$count++;
		print_it($id, $start, $stop, $scaf, $strand, $flag, $count);
	    }
	}
    }
}
#-----------------------------------------------------------------------------
sub add_lengths{

    my $stuff = shift;
    my $flag = shift;    
    my $id;
    my @col9 = split(/\;/, $stuff->[8]);
    
    if ( $flag eq 'i'){
        my @id_plus = split(/=/, $col9[1]);
	push(@{$ST_ST{$stuff->[0]}{$flag}{$id_plus[1]}{'stst'}},$stuff->[3]);
	push(@{$ST_ST{$stuff->[0]}{$flag}{$id_plus[1]}{'stst'}},$stuff->[4]);
	$ST_ST{$stuff->[0]}{'strand'}{$id_plus[1]} = $stuff->[6];
    }
    if ($flag eq 'p'){
	my @id_plus = split(/=/, $col9[1]);
	push(@{$P_GENES{$stuff->[0]}},$stuff->[3]);
	push(@{$P_GENES{$stuff->[0]}},$stuff->[4]);
    }
    else {
	my @id_plus = split(/=/, $col9[0]);
	push(@{$ST_ST{$stuff->[0]}{$flag}{$id_plus[1]}{'stst'}},$stuff->[3]);
	push(@{$ST_ST{$stuff->[0]}{$flag}{$id_plus[1]}{'stst'}},$stuff->[4]);
	$ST_ST{$stuff->[0]}{'strand'}{$id_plus[1]} = $stuff->[6];
    }   
}
#-----------------------------------------------------------------------------
sub parse_gff3{

    my $line = shift;
    my @stuff = split(/\t/, $line);

    if ($stuff[2] eq 'mRNA' && $opt_m){
	add_lengths(\@stuff, 'm');
    }
    if ($stuff[2] eq 'exon' && $opt_e){
	add_lengths(\@stuff, 'e');
    }
    if ($stuff[2] eq 'CDS' && $opt_c){
	add_lengths(\@stuff, 'c');
    }
    if ($stuff[2] eq 'five_prime_UTR' && $opt_u){
	add_lengths(\@stuff, 'f');
    }
    if ($stuff[2] eq 'three_prime_UTR' && $opt_u){
	add_lengths(\@stuff, 't');
    }
    if ($stuff[2] eq 'gene' && $opt_g){
	add_lengths(\@stuff, 'g');
    }
    if ($stuff[2] eq 'exon' && $opt_i){
	add_lengths(\@stuff, 'i');
    }
    if ($stuff[2] eq 'gene' && $opt_p){
	add_lengths(\@stuff, 'p');
    }
}

#-----------------------------------------------------------------------------
sub parse_fasta{

    my $line = shift;
    
    if ($line =~ /^>(\S+)/){
	$DEF = $1;
    }
    else {
	$SEQ{$DEF} .= $line;
    }
}
#-----------------------------------------------------------------------------
sub parse{

    my $file = shift;       
    my $flag = shift;
    
    my $fh = new FileHandle;
    $fh->open($file);
    
    my $i = 0;
    while (defined(my $line = <$fh>)){
	chomp($line);

	if ($flag eq 'g'){
	    last if $line =~ /^\#\#FASTA/;
	    next unless $line =~/^$opt_s/;
	    parse_gff3($line)  if $flag eq 'g';
	    
	}
	else {
	    parse_fasta($line);
	}
    }
    $fh->close();
}
#-----------------------------------------------------------------------------

