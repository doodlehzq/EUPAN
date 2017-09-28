#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-8-5
package adjAlign;
sub sam2bam{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t);
getopts("ht:");

my $usage="\nUsage: eupan sam2bam [options]  <mapping_directory> <output_directory> <samtools_directory>

eupan sam2bam is used to adjust mapping results including: 1)coverting sam to bam 2)sorting bam 3)merging 
bam 4)indexing bam.

The script will call samtools program, so the directory where samtools locates is needed. 

Necessary input description:
  mapping_directory      <string>     This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, One or more mapping results,
                                      *.sam, should exist.

  output_directory        <string>    Results will be output to this directory.To avoid 
                                      overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  QUAST_directory         <string>    samtools directory where executable samtools locates.

Options:
     -h                               Print this usage page.

     -t                   <int>       Threads used.
                                      Default: 1

";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$samtools_dir)=@ARGV;

#Detect executable quast.py
$samtools_dir.="/" unless($samtools_dir=~/\/$/);
my $samtools=$samtools_dir."samtools";
die("Error01: Cannot find samtools file in directory $samtools_dir\n") unless(-e $samtools);

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
}

#Read threads
my $thread_num=1;
$thread_num=$opt_t if(defined($opt_t));

#Adjust directory names and create output directory
$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

#Create output directory and sub-directories
mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);
my $err_data=$out_dir."data/";
mkdir($err_data);
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
#obtain *.sam files within the sample directory
    my @sam;
    opendir(ASS,$sd) || die("Error04: cannot open data directory: $sd\n");
    my @files=readdir(ASS);
    closedir ASS;
    foreach my $f (@files){
	if($f=~/(.+)\.sam$/){
	    push @sam,$1;
	}
    }
    if(@sam==0){
	print STDERR "Warnings: cannot find (*.bam) in $sd. Not processing!\n";
	next; 
    }

#create output directory for a sample
    my $sout=$out_data.$s."/";
    mkdir($sout) unless(-e $sout);
    my $serr=$err_data.$s."/";
    mkdir($serr) unless(-e $serr);
#generate command
    my $com="";
    foreach my $f (@sam){
#covert sam to bam
	$com.="$samtools view -\@ $thread_num -bS $sd$f.sam > $sout$f.bam\n";
#sort bam
	$com.="$samtools sort -\@ $thread_num -o $sout$f"."_sorted.bam $sout$f.bam\n";
    }
#merge
    if(@sam==1){
	foreach my $f (@sam){
	    $com.="mv ".$sout.$f."_sorted.bam $sout$s".".bam\n"
	}
    }
    else{
	$com.="$samtools merge -\@ $thread_num $sout$s".".bam";
	foreach my $f (@sam){
	    $com.=" ".$sout.$f."_sorted.bam";
	}
	$com.= "\n";;
    }
#index
    $com.="$samtools index $sout$s.bam\n";


#remove inter-files
    foreach my $f (@sam){
	my $rf="$sout$f.bam";
	$com.="rm $rf\n";
	$rf="$sout$f"."_sorted.bam";
	$com.="rm $rf\n";
    }
    system($com);
}
1;
}
1;
