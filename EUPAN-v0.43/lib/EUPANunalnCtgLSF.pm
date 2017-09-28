#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-11-5
package unalnCtg;

sub getUnaln{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_m);
getopts("hm");

my $usage="\nUsage: eupanLSF getUnalnCtg [options]  <assembly_directory> <QUAST_assess_directory> <output_directory>

eupanLSF getUnalnCtg is used to collect unaligned contigs.

Necessary input description:

  assembly_directory      <string>    This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, assembly results, including 
                                      file *_gc.contig, should exist.

  QUAST_assess_directory  <string>    This directory should contain many sub-directories 
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, quast assessment, including 
                                      directory file contigs_reports, should exist.

  output_directory        <string>    Results will be output to this directory.To avoid 
                                      overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

Options:
     -h                               Print this usage page.
     -m                               Use gap-closed contigs instead of raw contigs


";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($contig_dir,$quast_dir,$out_dir)=@ARGV;

#detect perl script to collect contigs for single sample
my $execp="getUnalnCtg";
my @path=split /:/,$ENV{PATH};
my $fpflag=0;
foreach my $p (@path){
    $p.="/".$execp;
    if(-e $p && -x $p){
	$fpflag=1;
	last;
    }
}
die("Executable getUnalnCtg cannot be found in your PATH!\n
") unless($fpflag);

#Adjust input directories
$contig_dir.="/" unless($contig_dir=~/\/$/);
$quast_dir.="/" unless($quast_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);
#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
}

my $thread_num=1;
my $suffix=".contig";
$suffix=".gcContig" if defined $opt_m;
#Create output directory and sub-directories
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

#Read samples
opendir(DATA,$contig_dir) || die("Error03: cannot open data directory: $contig_dir\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $cond=$contig_dir.$s."/";
    my $quad=$quast_dir.$s."/";
    unless(-d $cond && -d $quad){
	print STDERR "Warning: cannot find $s in $cond or $quad. Not processing.\n";
	next;
    }

#obtain *.contig file within the assembly directory
    my $contig_file="";
    opendir(ASS,$cond) || die("Error04: cannot open data directory: $cond\n");
    my @ass_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@ass_files){
	if($f=~/$suffix$/){
	    $contig_file=$f;
	    last;
	}
    }
    if($contig_file eq ""){
	print STDERR "Warnings: cannot find assembly file(*.contig or *.gcContig) in $cond\n";
	next; 
    }
    $contig_file=$cond.$contig_file;

#obtain *_gc.contig file within the assembly directory
    $quad.="contigs_reports/nucmer_output/";
    my $unalign_file="";
    opendir(ASS,$quad) || die("Error04: cannot open data directory: $quad\n");
    my @q_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@q_files){
	if($f=~/\.unaligned$/){
	    $unalign_file=$f;
	    last;
	}
    }
    if($unalign_file eq ""){
	print STDERR "Warnings: cannot find unaligned file(*.unaligned) in $quad\n";
	next; 
    }
    $unalign_file=$quad.$unalign_file;

#generate command
    my $outfile=$out_data.$s.".unaln.contig";
    my $com="$execp $contig_file $unalign_file $outfile $s";

#generate and submit job script
#************** Might be modified for different task submission system *******************
    my $job_file=$script_out."/".$s.".lsf";   #script_file
    my $err_file=$stderr_out."/".$s.".err";   #stderr_output_file
    my $out_file=$stdout_out."/".$s.".out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error05: Unable to create job file: $job_file\n");
    print JOB "\#BSUB -J $s","_un\n";              #job name
    #print JOB "\#BSUB -q fat\n";                     #queue name in the submission system
    print JOB "\#BSUB -o $out_file\n";               #stdout
    print JOB "\#BSUB -e $err_file\n";               #stderr
    print JOB "\#BSUB -n $thread_num\n";             #thread number
    print JOB "\#BSUB -R \"span[ptile=$thread_num]\"\n";
    print JOB "$com\n";                              #commands
    close JOB;
    system("bsub <$job_file");                       #submit job
#*****************************************************************************************
}
1;
}



sub mergeUnaln{
use strict;
use warnings;

my $usage="\nUsage: eupanLSF mergeUnalnCtg [options]  <unaln_directory> <output.fa>  

eupanLSF mergeUnalnCtg is used to merge unaligned contigs of each individuals to a single file.
";

die $usage if @ARGV!=2;
my ($dir,$out)=@ARGV;

die("$dir doesn't exist!
") unless -e $dir;
die("$dir is not a directory!
") unless -d $dir;
$dir.="/" unless $dir=~/\/$/;
opendir(DIR,$dir);
my @f=readdir(DIR);
closedir DIR;
my $com="cat";
foreach my $f (@f){
    next if $f=~/^\.+$/;
    my $nf=$dir.$f;
    $com.=" $nf";
}

$com.=" >$out";
system($com);
}
1;
