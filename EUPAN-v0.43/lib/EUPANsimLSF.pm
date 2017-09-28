#!/usr/bin/perl
#Created by Yue Zhao & Zhiqiang Hu, 2016-3-30
package sim;


sub simulation{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_n);
    getopts("hn:");
    my $usage="\nUsage: eupanLSF sim [options] <data_path> <output_directory>

sim is used to simulate the size of pan-genome and core-genome from gene presence-absence data.

Necessary input description:

  data_path   <string>           This path leads to gene.exist or geneFam.exist

  output_directory <string>      Both final output files and intermediate results 
                                 will be found in this directory. To avoid 
                                 overwriting of existing files. We kindly request
                                 that the output_directory should not exist. It is
                                 to say, this directory will be created by the 
                                 script itself.

Options:
     -h                          Print this usage page.
    
     -n            <int>         Specifies the number of random sampling times for
                                 simulation.
                                 default: 100
    
";

    die $usage if @ARGV!=2;
    die $usage if defined($opt_h);
    my ($gEfile,$out_dir)=@ARGV;

#Check existence of output directory
    if(-e $out_dir){
	die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
    }
    
#get simulation number
    my $number=100;
    if(defined($opt_n)){
	$number=$opt_n;
    }


#Adjust directory names and create output directory
    unless($out_dir=~/\/$/){
	$out_dir.="/";
    }
    mkdir($out_dir);

#read samples


##### read geneExist ########

open(EX,$gEfile);
chdir $out_dir;

open(OUT,">sim_out.txt");
my $i=0;
my @riceline;
my @geneName;
my @geneExt;
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
    my $u=0;
#    if(defined $knownCore{$item}){ 
#	foreach my $kk (@t){
#	    $kk=1;
#	}
#    }
   foreach my $kk (@t){
	$u+=$kk;
    }
    if($u==0){
	next;
    }
    push @geneName, $item;
    push @geneExt, \@t;
}
close EX;
###### simulation #######
print OUT "Time\tCore\tPan\n";
foreach (1..$number){
    my @r=frand(scalar(@riceline));
    my (@core,@pan);
    my ($c,$p)=(0,0);
#init
    for(my $i=0;$i<@geneExt;$i++){
	push @core,$geneExt[$i]->[$r[0]];
    }
    @pan=@core;
    for(my $i=0;$i<@pan;$i++){
	$c+=$core[$i];
    }
    $p=$c;
    print OUT "1\t$p\t$c\n";
    for(my $i=1;$i<@r;$i++){
	for(my $j=0;$j<@geneExt;$j++){
	    if($core[$j]==1 && $geneExt[$j]->[$r[$i]]!=1){
		$core[$j]=0;
		$c--;
	    }
	    if($pan[$j]!=1 && $geneExt[$j]->[$r[$i]]==1){
		$pan[$j]=1;
		$p++;
	    }
	}
	print OUT $i+1,"\t$c\t$p\n";
    }
}

sub frand{
    my @r;
    my %h;
    foreach (1..$_[0]){
	$h{$_-1}=1;
    }
    my @t=keys(%h);
    while(@t>0){
	my $i=int(rand(scalar(@t)));
	push @r,$t[$i];
	delete($h{$t[$i]});
	@t=keys(%h);
    }
    return @r;
}


close OUT;
    1;
}





sub pavPlot{
    use strict;
    use warnings;

    my $out_dir=$ARGV[1];
    unless($out_dir=~/\/$/){
	$out_dir.="/";
    }
	
    my $sim_input= "sim_out.txt";
    

#output qualities
    print STDERR "plotting simulation ...\n";
    plotSim($sim_input);


    sub plotSim{
        my ($file)=@_;
	my $exec="pav_plot.R";
	my @path=split /:/,$ENV{PATH};
	my $fpflag=0;
	foreach my $p (@path){
	    $p.="/".$exec;
	    if(-e $p && -x $p){
		$fpflag=1;
		$exec=$p;
		last;
	    }
       }
       if($fpflag){
          my $com="Rscript $exec $file 1>/dev/null 2>/dev/null\n";
          system($com);
          return 1;
       }
       else{
	   print STDERR ("pav_plot.R cannot be found in your PATH! Omit plotting!\n");
	   return 0;
       }
    }
    1;
}

1;

