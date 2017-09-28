#!/usr/bin/perl
use strict;
use warnings;

my $usage="$0 <fasta> <fasta.clstr> <conversion.txt> <out_fasta> <out.clstr>
";

die $usage if @ARGV!=5;

my ($fa,$facls,$con,$outfa,$outcls)=@ARGV;
my %h;
open(IN,$con);
while(<IN>){
    chomp;
    my ($a,$b)=split /\t/,$_;
    $h{$b}=$a;
}
close IN;

open(OUT,">$outfa");
open(IN,$fa);
while(<IN>){
    if(/^>([0-9]+)/){
	print OUT ">",$h{$1},"\n";
    }
    else{
	print OUT $_;
    }
}
close IN;
close OUT;

open(OUT,">$outcls");
open(IN,$facls);
while(<IN>){
    if(/^>/){
	print OUT $_;
    }
    else{
	my ($a,$b)=split /\.\.\./,$_;
	my($c,$d)=split /\>/,$a;
	print OUT $c,$h{$d},$b;
    }
}
close IN;


