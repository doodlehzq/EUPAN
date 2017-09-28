use strict;
use warnings;
package bam2cov;

sub bam2bed{
use strict;
use warnings;
my $usage="\nUsage: eupan bam2bed [options]  <bam_directory> <output_directory> 

This tool is used to calculate the covered region of the genome.
The outputs are covered fragments without overlap in 3-column .bed format. 

Necessary input description:

  bam_directory           <string>    This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, mapping result, a sorted .bam
                                      file, should exist.

  output_directory        <string>    Results will be output to this directory. To avoid 
                                      overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

";

die $usage if @ARGV<2;
my ($data_dir,$out_dir)=@ARGV;

#detect bam2cov

my $exec="bam2cov";
my @path=split /:/,$ENV{PATH};
my $fpflag=0;
foreach my $p (@path){
  $p.="/".$exec;
  if(-e $p && -x $p){
     $fpflag=1;
	last;
  }
}
die("Executable bam2cov cannot be found in your PATH!\n
") unless($fpflag);


#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
}

#Read threads
my $thread_num=1;

#Adjust directory names and create output directory
$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

#Create output directory and sub-directories
mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);

#Read samples
opendir(DATA,$data_dir) || die("Error03: cannot open data directory: $data_dir\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;

    my $sd=$data_dir.$s."/";
    next unless(-d $sd);
    print STDERR "Process sample $sd\n";
#obtain *.bam file within the sample directory
    my $bam_file=$sd.$s.".bam";
    unless(-e $bam_file){
	print STDERR "Warnings: cannot find bam file($bam_file) in $sd: skip this sample\n";
	next; 
    }

#create output directory for a sample
    my $sample_out.=$out_data."$s.bed";
#generate command
    my $com;
    $com="$exec $bam_file >$sample_out";
    system($com);
}
}
1;
