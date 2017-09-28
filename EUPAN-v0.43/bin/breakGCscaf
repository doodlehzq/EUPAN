#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-10-5
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_m);
getopts("m:");
my $usage="Usage: $0 [options] <gap-closed scaf> <gap-closed contig>
This script is to break gap-closed scaffolds into contigs at consecutive N.

Option:
   -m     <INT>    the number of consecutive N to be broken down
                   default: 10
";
die $usage if @ARGV!=2;

my $con_N=10;
$con_N=$opt_m if defined($opt_m);

break_scaffolds(@ARGV,$con_N);

sub break_scaffolds{
    my ($in,$out,$sepN)=@_;
    my $sep="";
    for(my $i=1;$i<$sepN;$i++){
	$sep.="N";
    }
    open(SCAF,$in)||die("Error09: cannot read file:$in\n");
    open(CONTIG,">$out")||die("Error10: cannot write file:$out\n");
    my $init=1;
    my ($header,$seq);
    while(my $line=<SCAF>){
	chomp $line;
	if($line=~/^>/){
	    if(!$init){
		split_and_print_seq($seq,$sep,$header);
		($header)=split /\s+/,$line;
		$seq="";
	    }
	    else{
		($header)=split /\s+/,$line;
		$seq="";
		$init=0;
	    }
	}
	else{
	    $seq.=$line;
	}
    }
    split_and_print_seq($seq,$sep,$header);
    close SCAF;
    close CONTIG;
}

sub split_and_print_seq{
    my ($seq,$sep,$header)=@_;
    my $i=0;
    my @tmp=split /\Q$sep\EN+/,$seq;
    foreach my $s (@tmp){
	next if length($s)==0;
	$i++;
	print CONTIG $header,"_",$i,"\n";
	for(my $j=0;$j<int(length($s)/100);$j++){
	    print CONTIG substr($s,100*$j,100),"\n";
	}
	print CONTIG substr($s,0-length($s)%100),"\n" if length($s)%100>0;
    }
}

