#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-7-2
package trim;
sub trimFastq{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_f $opt_t $opt_a $opt_s $opt_k $opt_i $opt_w $opt_n $opt_m $opt_j $opt_p $opt_q);
getopts("hf:t:a:s:k:i:w:n:m:j:p:q:");

my $usage="\nUsage: eupanLSF trim [options] <fastq_data_directory> <output_directory> <Trimmomatic_directory>

trim program is used to trim raw sequencing data to generate high-quality paired-end fastq data.

The script will call trimmomatic program and parameter files within trimmomatic
directory is also needed. So the directory where trimmomatic locates should be
given to the script as a necessary input.

Necessary input description:

  fastq_data_directory   <string>    This directory should contain many sub-directories
                                     named by sample names, such as CX101, B152,etc.
                                     In each sub-directory, there should be several 
                                     sequencing files ended by .fastq(or .fq) or .fastq.gz(or .fq.gz).

  output_directory       <string>    High-quality reads will be output to this directory. 
                                     To avoid overwriting of existing files. We kindly request
                                     that the output_directory should not exist. It is
                                     to say, this directory will be created by the 
                                     script itself.

  Trimmomatic_directory  <string>    Directory where trimmometic program locates.

Options:
     -h                              Print this usage page.

     -t                   <int>      thread number

     -a                   <string>   Adaptor file in fasta utilized by trimmomatic program.
                                     Default: trimmomatoc_dir/adapters/TruSeq3-PE-2.fa
     
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
     -p                   <33 or 64> Quality score version. 
                                     Default: 33 (phred+33)


     -i                   <int>      Parameter passed to Trimmomatic (LEADING or TRAILING). 
                                     Specifies the minimum  quality required to keep a base at the 
                                     head or tail of a read.
                                     Default: 20. 

     -w                   <int>      Parameter passed to Trimmomatic (SLIDINGWINDOW).specifies the 
                                     number of bases to average across.
                                     Default: 4.
     
     -n                   <int>      Parameter passed to Trimmomatic (quality for SLIDINGWINDOW).
                                     Default: 20.

     -m                   <int>      Parameter passed to Trimmomatic (MINLEN).Specifies the minimum 
                                     length of reads to be kept.
                                     Default: 35.

     -j                   <int>      Parameter passed to Trimmomatic (HEADCROP). The number of bases 
                                     to remove from the start of the read. If it is set to 0, no base
                                     will be removed from the head.
                                     Default: 0.
     -q                   <string>   The queue name for job submiting. 
                                     default: default queue
";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$trim_dir)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
}

#Detect executable fastqc
my $qc_exec;

#get thread number
my $thread_num=1;
if(defined($opt_t)){
    $thread_num=$opt_t;
}

unless($trim_dir=~/\/$/){
    $trim_dir.="/";
}

#get trimmomatic file
opendir(TRIM,$trim_dir)||die("Error: unable to open trimmomatic directory: $trim_dir\n");
my @trim_file=readdir(TRIM);
closedir TRIM;
my $trim_exec="";
foreach my $f (@trim_file){
    if($f=~/^trimmomatic.+\.jar$/){
	$trim_exec=$f;
	last;
    }
}
die("Error: unable to find \"trimmomatic*.jar\" in $trim_dir\n") if($trim_exec eq "");
$trim_exec=$trim_dir.$trim_exec;

#get trimmomatic adaptor file
my $trim_adaptor=$trim_dir."adapters/TruSeq3-PE-2.fa";
if(defined($opt_a)){
    $trim_adaptor=$opt_a;
}
die("Error: unable to find trimmomatic adaptor file: $trim_adaptor\n") unless(-e $trim_adaptor);

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

#read trimmomatic parameters
my $phred="phred33";
my $trim_leading=20;
my $trim_window=4;
my $trim_win_score=20;
my $trim_min=35;
my $trim_headcrop=0;
$trim_leading=$opt_i if(defined($opt_i));
$trim_window=$opt_w if(defined($opt_w));
$trim_win_score=$opt_n if(defined($opt_n));
$trim_min=$opt_m if(defined($opt_m));
$trim_headcrop=$opt_j if(defined($opt_j));
if(defined $opt_p){
    if($opt_p == "33"){
	next;
    }
    elsif($opt_p == "64"){
	$phred="phred64";
    }
    else{
	die "Unknown quality score version!
";
    }
}

#Adjust directory names and create output directory
unless($data_dir=~/\/$/){
    $data_dir.="/";
}
unless($out_dir=~/\/$/){
    $out_dir.="/";
}
mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);
#my $trim_log=$out_dir."trim_log/";
#mkdir($trim_log);
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
    my $sd=$data_dir.$s;
    if(-d $sd){
        #read sample directories
	opendir(RUN,$sd)|| die("Error: can not open directory: $sd\n");
	my @files=readdir(RUN);
	closedir RUN;
	my %fastq;
	foreach my $f (@files){
	    next if $f=~/^\./;
	    unless ($f=~/$suffix$/){
		print STDERR "Warning: file $f doesn't end with suffix: $suffix. This file won't be processed\n";
		next;
	    }
	    #put prefix of paired-end fastq files into %fastq
	    my @tmp=split /$linker/, $f;
	    pop @tmp;
	    my $nf=join "$linker",@tmp;
	    unless($f=~/^$nf\Q$linker\E[12]$suffix$/){
		die("Error: file $sd\/$f doesn't follow the \Q$linker\E[12]$suffix pattern\n");
	    }
	    $nf=$sd."/".$nf;
	    $fastq{$nf}=1 unless(defined($fastq{$nf}));
	}
        #generate commandline 
	my $tmp_out=$out_data.$s."/";
	mkdir($tmp_out);
#	my $log_out=$trim_log.$s."/";
#	mkdir($log_out);
	my $m=1;
	my $command="java -jar $trim_exec PE -$phred "; 
	foreach my $f (keys(%fastq)){
#	    my $c=$command."-trimlog $log_out".$m.".log ";  #logfile
	    my $c=$command.$f.$linker."1".$suffix." ".$f.$linker."2".$suffix." "; #input fastq
	    $c=$c.$tmp_out."paired_".$m."_1.fq ".$tmp_out."single_".$m."_1.fq ".$tmp_out."paired_".$m."_2.fq ".$tmp_out."single_".$m."_2.fq "; #output 4 files
	    $c=$c."ILLUMINACLIP:".$trim_adaptor.":2:30:10 "; #adaptor parameter file
	    if($trim_headcrop!=0){
		$c=$c."HEADCROP:".$trim_headcrop." ";  #headcrop
	    }
	    $c=$c."LEADING:".$trim_leading." TRAILING:".$trim_leading." "; #leading and trailing
	    $c=$c."SLIDINGWINDOW:".$trim_window.":".$trim_win_score." ";   #sliding window
	    $c=$c."MINLEN:".$trim_min;
	    # gzip fastq files
	    my $gcom1="gzip $tmp_out"."paired_".$m."_1.fq";
	    my $gcom2="gzip $tmp_out"."paired_".$m."_2.fq";
	    my $gcom3="gzip $tmp_out"."single_".$m."_1.fq";
	    my $gcom4="gzip $tmp_out"."single_".$m."_2.fq";

#************** Might be modified for different task submission system *******************
	    my $job_file=$script_out."/".$s.".".$m.".lsf";   #script_file
	    my $err_file=$stderr_out."/".$s.".".$m.".err";   #stderr_output_file
	    my $out_file=$stdout_out."/".$s.".".$m.".out";   #stdout_output_file
	    #create job script
	    open(JOB,">$job_file")||die("Error: Unable to create job file: $job_file\n");
	    print JOB "\#BSUB -J $s\n";                      #job name
	    print JOB "\#BSUB -q $opt_q\n" if defined $opt_q;   #queue name in the submission system
	    print JOB "\#BSUB -o $out_file\n";               #stdout
	    print JOB "\#BSUB -e $err_file\n";               #stderr
	    print JOB "\#BSUB -n $thread_num\n";             #thread number
	    print JOB "\#BSUB -R \"span[ptile=$thread_num]\"\n"; #use the same node
	    print JOB "$c\n";                                #commands
	    print JOB "$gcom1\n";                            #commands
	    print JOB "$gcom2\n";                            #commands
	    print JOB "$gcom3\n";                            #commands
	    print JOB "$gcom4\n";                            #commands
	    close JOB;
	    print STDERR "submit jobs for run: $f\n";
	    system("bsub <$job_file");                       #submit job
#*****************************************************************************************
	    $m++;
	}
    }
    else{
	print STDERR "Warning: $sd is not a directory! =>  Not processed!\n";
    }
}
1;
}
1;
