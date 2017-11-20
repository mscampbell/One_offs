#!/usr/bin/perl -w 
use strict;
use lib ('/home/mcampbell/lib');

my $FILE = $ARGV[0];
my $outfile = $FILE;
$outfile =~ s/gff$/gtf/;
my $exe = '/home/mcampbell/applications/cufflinks-2.2.1.Linux_x86_64/gffread';
my $options = '-T -g /home/mcampbell/projects/EVM_MAKER/bench_marking/drosophila/assembly/dmel_3R.fasta';

my $command = $exe;
$command .= " $FILE";
$command .= " $options";
$command .= " -o $outfile";

system($command);
