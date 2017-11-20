#!/usr/bin/perl -w 
use strict;
use Getopt::Std;
use vars qw($opt_d);
getopts('d');

use warnings;
use lib ('/Users/mcampbell/mcampbell/lib');
use PostData;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "\n\nscript to pull out the maker genes that are not 
in the intersection of maker and ensembl genes and see if they
have pfam domains \n\n             
OPTIONS: -d <prints the top 100 pfam ids for the Maker only genes>\n
USAGE: lamprey_script.pl <file1.rbh> <file2.maker.gff>\n\n";

die($usage) unless $ARGV[1];

my $FILE1 = $ARGV[0];
my $FILE2 = $ARGV[1];
my %GFF3;
my $RBHPFAM=0;
my %INDEX;
my $MEINTER=0;
my $TPFAM=0;
my $TMG=0;
my $MAKER_only=0;
my @MAKER_ONLY;
my %PFDS;
my %PFDS2;
my @RBHGENES;

parse_rbh($FILE1);
parse_gff3($FILE2);
report();
report2();
pfd() if $opt_d;
#print "Maker only\n";
#pfds();
#print "\nRBH\n";
#pfds2();
print "\n\nThank you for choosing Mike\'s lamprey_script.pl\n\n";

#-----------------------------------------------------------------------------
#---------------------------------- SUBS -------------------------------------
#-----------------------------------------------------------------------------
#prints the top 100 pfam ids for the Maker only genes
sub pfds{
    my $number = 0;
    foreach my $dom (sort {$PFDS{$b} <=> $PFDS{$a}} keys %PFDS){
        print $dom."\t".$PFDS{$dom}."\n" unless $number >= 100;
        $number++;
    }
}
#-----------------------------------------------------------------------------
#prints the top 100 pfam ids for the RBH genes                                                   
sub pfds2{
    my $number = 0;
    foreach my $dom (sort {$PFDS2{$b} <=> $PFDS2{$a}} keys %PFDS2){
        print $dom."\t".$PFDS2{$dom}."\n" unless $number >= 100;
        $number++;
    }
}

#-----------------------------------------------------------------------------
sub parse_rbh{
    my $file = shift;
    
    open (my $FH, '<', $file) || die "Cant\'t find $file\n";
    
    while (defined(my $line = <$FH>)){
        chomp($line);
	
	my ($ens, $maker) = split(/\s+/,$line);
	$INDEX{$maker}=$ens;
    }
    close($FH);
}
#-----------------------------------------------------------------------------
sub parse_gff3{
    my $file = shift;
    
    open (my $FH, '<', $file) || die "Cant\'t find $file\n";
    
    while (defined(my $line = <$FH>)){
        chomp($line);
	last if $line =~ /^\#\#FASTA/;
	next if $line =~ /^\#/;
	next if $line =~ /^\s*$/;    
	next if $line =~ /^\>/;    
	my @GFF_col = split(/\t/, $line);
	
	if ($GFF_col[2] eq 'mRNA'){
	    my @mrna9 = split(/\;/, $GFF_col[8]); 
	    $mrna9[0] =~ s/ID=//;
	    $mrna9[0] =~ s/\-RA//;
	    foreach my $pair (@mrna9){
		my ($k, $v) = split(/\=/,$pair);
		$GFF3{$mrna9[0]}{$k}=$v;
		$GFF3{$mrna9[0]}{ensemble_hit}=0;
	    }           
	}
	
    }
    close($FH);
}
######################### example data structure ########################
#{PMZ_0024127-RA} = HASH
#{Alias} = maker-scaffold_23256.1-1073-augustus-gene-0.0-mRNA-1
#{Dbxref} = InterPro:IPR001593,InterPro:IPR018281,PANTHER:PTHR11830,Pfam:PF01015,Prosite:PS01191
#{Name} = PMZ_0024127-RA
#{Note} = Similar to rps3a: 40S ribosomal protein S3a (Danio rerio)
#{Ontology_term} = GO:0003735,GO:0005622,GO:0005840,GO:0006412
#{PMZ_0024127-RA} = UNDEFINED VALUE
#{Parent} = PMZ_0024127
#{_AED} = 0.00
#{_QI} = 0|-1|0|1|-1|1|1|0|274

#-----------------------------------------------------------------------------
sub pfd{
    foreach my $x (@MAKER_ONLY){
	if (defined($GFF3{$x}{Dbxref})){
	    my @stuff = split(/\,/,$GFF3{$x}{Dbxref});
	    foreach my $y (@stuff){
		if ($y =~ /^Pfam:/){
		    my ($k, $v) = split(/\:/,$y);
		    $PFDS{$v}++; 
		}
	    }    
	}
    }
    foreach my $x (@RBHGENES){
        if (defined($GFF3{$x}{Dbxref})){
            my @stuff = split(/\,/,$GFF3{$x}{Dbxref});
            foreach my $y (@stuff){
                if ($y =~ /^Pfam:/){
                    my ($k, $v) = split(/\:/,$y);
                    $PFDS2{$v}++;
                }
            }
        }
    }


}
#-----------------------------------------------------------------------------
sub report{
    my @ikeys= keys %INDEX;
    foreach my $key (@ikeys){
	$GFF3{$key}{ensemble_hit}++;
	$MEINTER++;
	my $x = $GFF3{$key}{ensemble_hit};
	my $y = $GFF3{$key}{Dbxref};
	if ($x == 1 && $y){
	    push(@RBHGENES, $key);
	    $RBHPFAM++;
# prints the list of RBH hits with pfams
#	    print "EandM ".$key . "\t". $y . "\n";
	}	
    }
    
    print "\n\nMAKER Ensemble intersection:\t". $MEINTER . "\n";
    foreach my $key (keys %GFF3){
	if ($GFF3{$key}{Dbxref}){
	    $TPFAM++
	    }
	
    } 
    print "total maker genes with Pfam domains:\t". $TPFAM . "\n";
    foreach my $key (keys %GFF3){
	$TMG++
	}
    
    print "total maker genes:\t". $TMG . "\n"
    }

#-----------------------------------------------------------------------------
sub report2{
    foreach my $mkey (keys %GFF3){                                                       
	my $g = $GFF3{$mkey}{ensemble_hit};                                              
	my $f = $GFF3{$mkey}{Dbxref};                                                    
	if ($g == 0 && $f){
	    push(@MAKER_ONLY, $mkey);                                                              
	    $MAKER_only++;	 
# prints the maker only genes with pfams
#	    print "M_only ".$mkey . "\t". $f . "\n";
	}
	
    }
    print "total maker only genes with pfam domains:\t". $MAKER_only . "\n";
    print "RBH genes with pfam domains:\t". $RBHPFAM . "\n\n";

}
#-----------------------------------------------------------------------------
