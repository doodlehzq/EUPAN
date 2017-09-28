#!/usr/bin/perl
package gExist;

sub checkGeneExist{
    use strict;
    use warnings;

my $usage="eupan geneExist gene_file cds_file min_gene_cov min_cds_cov >output

Inputs:

gene_file     <string>   gene body coverage file
cds_file      <string>   CDS coverage file
min_gene_cov  <float>    minimum gene body coverage
min_cds_cov   <float>    minimum CDS coverage

";
die $usage if @ARGV!=4;
my ($gf,$cf,$gcov,$ccov)=@ARGV;

open(GF,$gf);
open(CF,$cf);
my $i=0;
while(<GF>){
    chomp;
    my $c=<CF>;
    chomp $c;
    $i++;
    if($i==1){
	print $_,"\n";
	next;
    }
    else{
	my @t1=split /\t/,$_;
	my @t2=split /\t/,$c;
	print $t1[0];
	for(my $i=1;$i<@t1;$i++){
	    if($t1[$i]>=$gcov && $t2[$i]>=$ccov){
		print "\t1";
	    }
	    else{
		print "\t0";
	    }
	}
	print "\n";
    }
}
close GF;
close CF;
}

sub subsetSample{
use strict;
use warnings;

my $usage="eupan subSample gene_existence_matrix sample_list >output
";

die $usage if @ARGV!=2;

my ($matrix,$rice)=@ARGV;

my %h;
open(IN,$rice);
while(<IN>){ 
    chomp;
    $h{$_}=1;
}
close IN;
my $i=0;
my @n;

open(IN,$matrix);
while(<IN>){
    $i++;
    my @t=split /\t/,$_;
    if($i==1){
	foreach (@t){
	    if(defined $h{$_}){
		push @n,1;
	    }
	    else{
		push @n,0;
	    }
	    $n[0]=1;
	}
	output_line(\@n,\@t);
    }
    else{
	output_line(\@n,\@t);
    }
}
close IN;


sub output_line{
    my ($nad,$ad)=@_;
    print $ad->[0];
    for(my $i=1;$i<@{$ad};$i++){
	print "\t",$ad->[$i] if $nad->[$i];
    }
    print "\n";
}
}

sub gE2gfE{
use strict;
use warnings;

my $usage="eupan gFamExist <geneExist.txt> <geneFam.info> > geneFamExist.txt
";

die $usage if @ARGV!=2;

my ($gEfile,$finfo)=@ARGV;

########## read group #######
my @group;
my @gr_index;

open(GF,$finfo);
while(<GF>){
    chomp;
    my @t=split /\s+/,$_;
    my $gid=shift @t;
    $gid=substr($gid,0,length($gid)-1);
    push @group,$gid;
    push @gr_index,\@t;
}
close GF;
############ end ############

##### read geneExist ########

open(EX,$gEfile);
my $i=0;
my @riceline;
my %geneExt;
while(<EX>){
    $i++;
    chomp;
    if($i==1){
	@riceline=split /\t/,$_;
	shift @riceline;
	next;
    }
    my @t=split /\t/,$_;
    my $item=shift @t;
    $geneExt{$item}=\@t;
}
close EX;
print STDERR "Gene family number:",scalar(@group),"\n";
print STDERR "Gene number:",scalar(keys(%geneExt)),"\n";
print STDERR "Rice line number:",scalar(@riceline),"\n";
############ end ############

###### gene fam exist #######

print "Gene\t",join("\t",@riceline),"\n";

for($i=0;$i<@group;$i++){
    print $group[$i];
    for(my $j=0;$j<@riceline;$j++){
	my $f=0;
	foreach my $k (@{$gr_index[$i]}){
	    if(!defined $geneExt{$k}){
		print STDERR $k,"\n";
		die;
	    }
	    if(!defined $geneExt{$k}->[$j]){
		print STDERR "i",$i,"\t",scalar(@{$gr_index[$i]}),"j$j","\n";
		next;
	    }
	    if($geneExt{$k}->[$j]==1){
		$f=1;
		last;
	    }
	}
	if($f==1){
	    print "\t1";
	}
	else{
	    print "\t0";
	}
    }
    print "\n";
}
}

1;
