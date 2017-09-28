use strict;
use warnings;
package bam2cov;
use Getopt::Std;
use vars qw($opt_q);
getopts("q:");

sub bam2bed{
use strict;
use warnings;
my $usage="\nUsage: eupanSLURM bam2bed [options]  <bam_directory> <output_directory> 

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

     -q            <string>      The queue name for job submiting. 
                                 default: default queue
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
#************** Might be modified for different task submission system *******************
my $job_out=$out_dir."job";       
mkdir($job_out);
my $script_out=$job_out."/scripts"; #job script directory
mkdir($script_out);
my $stderr_out=$job_out."/err";     #stdout directory
mkdir($stderr_out);
my $stdout_out=$job_out."/out";     #sdterr directory
mkdir($stdout_out);
#*****************************************************************************************

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
    $com="$exec $bam_file >$sample_out\n";

#generate and submit job script
#************** Might be modified for different task submission system *******************
    my $job_file=$script_out."/".$s."_bam2bed.slurm";   #script_file
    my $err_file=$stderr_out."/".$s."_bam2bed.err";   #stderr_output_file
    my $out_file=$stdout_out."/".$s."_bam2bed.out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error05: Unable to create job file: $job_file\n");
	    print JOB "\#!/bin/bash\n";
    print JOB "\#SBATCH --job-name=$s","_bam2bed\n";              #job name
	    print JOB "\#SBATCH -p $opt_q\n" if defined $opt_q;   #queue name in the submission system
    print JOB "\#SBATCH --output=$out_file\n";               #stdout
    print JOB "\#SBATCH --error=$err_file\n";               #stderr
    print JOB "\#SBATCH -n $thread_num\n";             #thread number
    print JOB "$com\n";                              #commands
    close JOB;
    system("sbatch $job_file");                       #submit job
#*****************************************************************************************


}
}
1;
