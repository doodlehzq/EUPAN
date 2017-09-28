#!/usr/bin/perl
package pTpG;
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_e);
getopts("e");

sub pTpG{
my $usage = "eupanLSF perTranPerGene <input_file.gtf> <output>

This tool is to obtain the longest trancript of each gene.

Option:
   -e    Check \"exon\" length instead of check \"CDS\" length
         Note \"exon\" or \"CDS\" should exist in the 3rd column of the input file

";
die $usage if @ARGV!=2;
my ($file,$out) = @ARGV;
my ($Vec, $num) = @{Read_GTF($file)};

open(OUT,">$out");
foreach my $d (@{$Vec}){
    my $target=$d->[0];
    my $len=getLen($d->[0]);
    for(my $i=1;$i<@{$d};$i++){
	my $nlen=getLen($d->[$i]);
	if($nlen>$len){
	    $len=$nlen;
	    $target=$d->[$i];
	}
    }
    if($len!=0){
	OutputTran($target);
    }
}

exit;


########## sub routines ###########
sub getLen{
    my ($d)=@_;
    my $len=0;
    foreach my $k (@{$d}){
	if(defined $opt_e){
	    if($k->{type} eq "exon"){
		$len+=$k->{end}-$k->{start}+1;
	    }
	}
	else{
	    if($k->{type} eq "CDS" ||$k->{type} eq "stop_codon"){
		$len+=$k->{end}-$k->{start}+1;
	    }
	}
    }
    return $len;
}

sub OutputTran{              #\@tvec,$name ,output file
    my ($ad) = @_;
    my $d1;
    for($d1=0;$d1<@{$ad};$d1++){
	print OUT "${${$ad}[$d1]}{chr}\t";
	print OUT "${${$ad}[$d1]}{source}\t";
	print OUT "${${$ad}[$d1]}{type}\t";
	print OUT "${${$ad}[$d1]}{start}\t";
	print OUT "${${$ad}[$d1]}{end}\t";
	print OUT "${${$ad}[$d1]}{score}\t";
	print OUT "${${$ad}[$d1]}{sym}\t";
	print OUT "${${$ad}[$d1]}{phase}\t";
	print OUT "${${$ad}[$d1]}{gid}; ";
	print OUT "${${$ad}[$d1]}{tid};\n";
    }
}

sub Read_GTF{
    my ($file) = @_;
    my @gvec;
    my @re;
    my $begin = 1;
    my $pre_tid;
    my $pre_gid;
    my $i = 0;
    my $j = 0;
    my $k = 0;
    my @temp;
    open(IN,$file)||die "open $file error:$!\n";
    while(<IN>){
	chomp;
	if(/^[0-9a-zA-Z]+/){

	    if($begin == 1){
		@temp = split("\t",$_);
		$temp[8]=~/(gene_id \"[^\"]+\"); (transcript_id \"[^\"]+\")/;
		$pre_tid = $2;
		$pre_gid = $1;

		$gvec[$i] = init_gene();
		${$gvec[$i]}[$j] = init_tran();

		${${$gvec[$i]}[$j]}[$k] =
			  new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
		$begin = 0;
	    }
	    else{
		@temp = split("\t",$_);
		$temp[8] =~ /(gene_id \"[^\"]+\"); (transcript_id \"[^\"]+\")/;

		if($1 eq $pre_gid){
		    if($2 eq $pre_tid){
			$k++;
			${${$gvec[$i]}[$j]}[$k] =
				  new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
		    }
		    else{
			$j++;
			${$gvec[$i]}[$j] = init_tran();
			$k = 0;
			${${$gvec[$i]}[$j]}[$k] = 
				  new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
			$pre_tid = $2;
		    }
		}
		else{
		    $i++;
		    $gvec[$i] = init_gene();
		    $j = 0;
		    ${$gvec[$i]}[$j] = init_tran();
		    $k = 0;
		    ${${$gvec[$i]}[$j]}[$k] = 
			      new_sf($temp[0],$temp[1],$temp[2],$temp[3],$temp[4],$temp[5],$temp[6],$temp[7],$1,$2);
		    $pre_gid =$1;
		    $pre_tid = $2;
		}
	    }
	}
    }
    close IN;
    $re[0] = \@gvec;
    $re[1] = scalar(@gvec);
    return \@re;
}

sub new_sf{             
    my ($chr,$src,$type,$start,$end,$score,$sym,$phase,$gid,$tid)=@_;
    my %hash = (
	"chr" => $chr,
	"source" => $src,
	"type" => $type,
	"score" => $score,
	"gid" => $gid,
	"tid" => $tid,
	"start" => $start,
	"end" => $end,
	"sym" => $sym,
	"phase" => $phase,
    );
    return \%hash;
}

sub init_gene{
    my @a;
    return \@a;
}
sub init_tran{
    my @a;
    return \@a;
}
1;
}
1;
