#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-8-5
package align;
sub map{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_s $opt_n $opt_m $opt_f $opt_k $opt_q);
getopts("hs:n:m:t:f:k:q:");

my $usage="\nUsage: eupanSLURM align [options]  <fastq_data_directory> <output_directory> <Mapping_tool_directory> <alignment_index>

eupanSLURM align is used to map high-quality reads to a reference on large scale.

The script will call mapping program (bwa mem or bowtie2), so the directory where mapping tool locates is needed. 

Necessary input description:

  fastq_data_directory    <string>    This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, there should be several 
                                      sequencing files ended by .fq(.gz) or .fastq(.gz).

  output_directory        <string>    Alignment results will be output to this directory.
                                      To avoid overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  Mapping_tool_directory  <string>    Directory where bwa or bowtie2 program locates.

  Alignment_index         <string>    bowtie2 index built from bowtie2-build program;
                                      or bwa index built from bwa index.   

Options:
     -h                              Print this usage page.

     -f                   <string>   Select a mapping tool. Can be \"bwa\" or \"bowtie2\".
                                     Default: bwa

     -t                   <int>      Threads used.
                                     Default: 1

     -s                   <string>   Suffix of the fastq_file. Check your sequencing data and
                                     change it if needed.
                                     Default: \".fq.gz\"

     -k                   <string>   Linker for paired_end identifer. Paired-end fastq file
                                     should end with *1suffix or *2suffix, where suffix is
                                     \".fq.gz\"( or \".fastq\", etc. See -s option) and * is the
                                     linker such as \"_\".As an example, the file should 
                                     be like CX123_1.fq.gz (linker is \"_\", suffix is \".fq.gz\")
                                     or BX125_R1.fastq(linker is \"_R\", suffix is \".fastq\")
                                     Default: \"_\"

     -m                   <int>      min insertion length, for bowtie2 only 
                                     Default: 0
   
     -n                   <int>      max insertion length, for bowtie2 only
                                     Default: 1000

     -q            <string>      The queue name for job submiting. 
                                  default: default queue

";

die $usage if @ARGV!=4;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$map_dir,$map_index)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.");
}

#select tools
my $tool="bwa";
if(defined($opt_f)){
    $tool=$opt_f;
}
die("Unknown tool: $tool\nPlease use -f bowtie2 or -f bwa\n") unless($tool eq "bwa" || $tool eq "bowtie2");

#Detect executable bowtie2 or bwa
$map_dir.="/" unless($map_dir=~/\/$/);
my $exec;
if($tool eq "bwa"){
    $exec=$map_dir."bwa";
    die("Cannot found executable bwa file in directory $map_dir\n") unless(-e $exec);
}
else{
    $exec=$map_dir."bowtie2";
    die("Cannot found executable bowtie2 file in directory $map_dir\n") unless(-e $exec);
}

#read threads
my $thread_num=1;
if(defined($opt_t)){
    $thread_num=$opt_t;
}

#read fastq suffix
my $suffix=".fq.gz";
if(defined($opt_s)){
    $suffix=$opt_s;
}

#read linker
my $linker="_";
if(defined($opt_k)){
    $linker=$opt_k;
}

#read  parameters
my ($min,$max)=(0,1000);
if(defined($opt_m)){
    $min=$opt_m;
}
if(defined($opt_n)){
    $max=$opt_n;
}



#Adjust directory names and create output directory

$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);


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


#read samples
opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample

foreach my $s (@sample){
    next if $s=~/^\./;
    my $sd=$data_dir.$s."/";
    print STDERR "Warning: $sd isn't a directory! => Not processed.\n" unless -d $sd;
    next unless -d $sd;
    #read sample directories
    opendir(RUN,$sd)|| die("Error: can not open directory: $sd\n");
    my @files=readdir(RUN);
    close RUN;
    my %fastq;
    foreach my $f (@files){
	next if $f=~/^\./;
	next if $f=~/^single/;
	unless ($f=~/$suffix$/){
	    print STDERR "Warning: file $f doesn't end with suffix: $suffix. This file won't be processed";
	}
	#put prefix of paired-end fastq files into %fastq
	my @tmp=split /$linker/, $f;
	pop @tmp;
	my $nf=join "$linker",@tmp;
	unless($f=~/^$nf\Q$linker\E[12]$suffix$/){
	    print STDERR "Warnings: file $sd\/$f doesn't follow the \Q$linker\E[12]$suffix pattern! =>\n";
	}
	$nf=$sd.$nf;
	$fastq{$nf}=1 unless(defined($fastq{$nf}));
    }

    my $tmp_out=$out_data.$s."/";
    mkdir($tmp_out) unless(-e $tmp_out);
    my $irr=0;
    foreach my $k (keys(%fastq)){
	$irr++;
	my $forward=$k.$linker."1".$suffix;
	my $reverse=$k.$linker."2".$suffix;
	my $com;
	if($tool eq "bowtie2"){
	    $com=$exec." -x $map_index"." -1 $forward -2 $reverse -S $tmp_out".$irr."_paired.sam"." -I $min -X $max -p $thread_num";
	}
	else{
	    $com="$exec mem -t $thread_num $map_index $forward $reverse > $tmp_out".$irr."_paired.sam";
	}
#************** Might be modified for different task submission system *******************
	my $job_file=$script_out."/".$s."_$irr.slurm";   #script_file
	my $err_file=$stderr_out."/".$s."_$irr.err";   #stderr_output_file
	my $out_file=$stdout_out."/".$s."_$irr.out";   #stdout_output_file
	#create job script
	open(JOB,">$job_file")||die("Error: Unable to create job file: $job_file\n");
	    print JOB "\#!/bin/bash\n";
	print JOB "\#SBATCH --job-name=$s","_$irr"."_map\n";              #job name
	print JOB "\#SBATCH -p $opt_q\n" if defined $opt_q;   #queue name in the submission system
	print JOB "\#SBATCH --output=$out_file\n";               #stdout
	print JOB "\#SBATCH --error=$err_file\n";               #stderr
#	print JOB "\#SBATCH -p fat\n";               #stderr
	print JOB "\#SBATCH -n $thread_num\n";             #thread number
	print JOB "\#SBATCH --ntasks-per-node=$thread_num\n";
	print JOB "$com\n";                              #commands
	close JOB;
	system("sbatch $job_file");                       #submit job
#*****************************************************************************************

    }
}
1;
}
1;
