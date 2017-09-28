#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-7-4
package qualStaSLURM;
sub checkQual{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_f $opt_t $opt_q);
    getopts("hf:t:q:");

    my $usage="\nUsage: eupanSLURM qualSta [options] <data_directory> <output_directory>

qualSta is used to check qualities of .fastq/.fastq.gz files on a large scale.

The script will call fastqc program, so please make sure fastqc is in your
PATH, or you need to use -f option to tell the script where it locates.

Necessary input description:

  data_directory   <string>      This directory should contain many sub-directories
                                 named by sample names, such as CX101, B152,etc.
                                 In each sub-directory, there should be several 
                                 sequencing files ended by .fastq or .fastq.gz.

  output_directory <string>      Both final output files and intermediate results 
                                 will be found in this directory. To avoid 
                                 overwriting of existing files. We kindly request
                                 that the output_directory should not exist. It is
                                 to say, this directory will be created by the 
                                 script itself.

Options:
     -h                          Print this usage page.

     -f            <string>      The directory where the executable file (fastqc) 
                                 locate. If this option isn't given, we assume 
                                 that it is in your PATH.
     
     -t            <int>         Specifies the number of files which can be processed
                                 simultaneously. This parameter is sent to fastqc 
                                 program. It is recommended to set as the number of 
                                 files within each sample. Pay attention that the
                                 machine should have this number of threads.
                                 default: 1 

     -q            <string>      The queue name for job submiting. 
                                 default: default queue
";

    die $usage if @ARGV!=2;
    die $usage if defined($opt_h);
    my ($data_dir,$out_dir)=@ARGV;
    
#Check existence of output directory
    if(-e $out_dir){
	die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.");
}
    
#Detect executable fastqc
    my $qc_exec="fastqc";
    if(defined($opt_f)){
	unless($opt_f=~/\/$/){
	    $opt_f.="/";
	}
	$qc_exec=$opt_f."fastqc";
	unless(-e $qc_exec && -x $qc_exec){
	    die("Error: $qc_exec doesn't exist or it is not executable!\n");
	}
    }
    else{
	my @path=split /:/,$ENV{PATH};
	my $fpflag=0;
	foreach my $p (@path){
	    $p.="/fastqc";
	    if(-e $p && -x $p){
		$fpflag=1;
		last;
	    }
	}
	unless($fpflag){
	    die("Error: fastqc doesn't exist in PATH!\nTry -f option!\n");
	}
    }

#get thread number
    my $thread_num=1;
    if(defined($opt_t)){
	$thread_num=$opt_t;
    }


#Adjust directory names and create output directory
    unless($data_dir=~/\/$/){
	$data_dir.="/";
    }

    unless($out_dir=~/\/$/){
	$out_dir.="/";
    }
    mkdir($out_dir);
    
#read samples
    opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
    my @sample=readdir(DATA);
    closedir DATA;

#create directory for fastqc result
    my $fastqc_out=$out_dir."fastqc_output";
    mkdir($fastqc_out);

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
	my $sd=$data_dir.$s;
	if(-d $sd){
        #read sample directories
	    opendir(RUN,$sd)|| die("Error: can not open directory: $sd\n");
	    my @files=readdir(RUN);
	    close RUN;
	    my @fastq;
	    foreach my $f (@files){
		next if $f=~/^\./;
		$f=$sd."/".$f;
		if($f=~/\.fastq$/ || $f=~/\.fastq\.gz$/ || $f=~/\.fq$/ || $f=~/\.fq\.gz$/){
                #put fastq files into @fastq
		    push @fastq, $f;
		}
		else{
		    print STDERR "Warning: $f is not a .fastq or .fastq.gz file! =>  Not processed!\n";
		}
	    }
        #generate commandline
	    my $command=$qc_exec;
	    foreach my $f (@fastq){
		$command.= " $f";
	    }
	    my $tmp_out=$fastqc_out."/".$s;
	    mkdir($tmp_out);
	    $command.= " -o $tmp_out -t $thread_num -q --extract";

#************** Might be modified for different task submission system *******************
	    my $job_file=$script_out."/".$s.".slurm";   #script_file
	    my $err_file=$stderr_out."/".$s.".err";   #stderr_output_file
	    my $out_file=$stdout_out."/".$s.".out";   #stdout_output_file
	    #create job script
	    open(JOB,">$job_file")||die("Error: Unable to create job file: $job_file\n");
	    print JOB "\#!/bin/bash\n";
	    print JOB "\#SBATCH --job-name=$s\n";              #job name
	    print JOB "\#SBATCH -p $opt_q\n" if defined $opt_q;   #queue name in the submission system
	    print JOB "\#SBATCH --output=$out_file\n";       #stdout
	    print JOB "\#SBATCH --error=$err_file\n";       #stderr
	    print JOB "\#SBATCH -n $thread_num\n";     #thread number
	    print JOB "$command\n";                  #commands
	    close JOB;
	    print STDERR "submit jobs for sample: $s\n";
	    system("sbatch $job_file");               #submit job
#*****************************************************************************************
	}
	else{
	    print STDERR "Warning: $sd is not a directory! =>  Not processed!\n";
	}
    }
    1;
}
sub mergeFastqc{
    use strict;
    use warnings;
    die "Usage: eupanSLURM mergeQualSta <directory>
directory  <string>  The exact output of \"eupanSLURM qualSta\".
" if @ARGV<1;
    my $out_dir=$ARGV[0];
    unless($out_dir=~/\/$/){
	$out_dir.="/";
    }
    my $fastqc_out=$out_dir."fastqc_output";
    
    opendir(SAMPLE,$fastqc_out) || die("Error: Unable to find fastqc directory: $fastqc_out\n");
    my @sample=readdir(SAMPLE);
    closedir SAMPLE;


#Collect fatsqc results
    my @quality_file;
    my @quality_sample;
#headers
    my @header=qw(Basic_Statistics Per_base_sequence_quality Per_tile_sequence_quality Per_sequence_quality_scores
 Per_base_sequence_content Per_sequence_GC_content Per_base_N_content
 Sequence_Length_Distribution Sequence_Duplication_Levels Overrepresented_sequences Adapter_Content Kmer_Content);

    foreach my $s (@sample){
	next if $s=~/^\./;
	my $sd=$fastqc_out."/".$s;
	opendir(FASTQC,$sd)||die("Unable to open fastqc output directory for sample $s: $sd\n");
	my @items=readdir(FASTQC);
	closedir FASTQC;
	my @sample_info_tmp;
	foreach my $f (@items){
	    next if $f=~/^\./;
	    next if $f=~/\.zip$/;
	    next if $f=~/\.html$/;
	    my $fd=$sd."/".$f."/"."summary.txt";
	    open(FILE,$fd) ||die("Unable to open fastqc output file: $fd\n");
	    my @tmp=<FILE>;
	    close FILE;
	    my @item_info;
	    push @item_info,$s."_".$f;
	    foreach my $t (@tmp){
		my ($r)=split /\s+/,$t;
		if($r eq "PASS"){
		    push @item_info,1;
		}
		elsif($r eq "WARN"){
		    push @item_info,0;
		}
		elsif($r eq "FAIL"){
		    push @item_info,-1;
		}
		else{
		    die("Error: unknown item in file $fd: $r\n");
		}
	    }
	    push @sample_info_tmp,\@item_info;
	    push @quality_file,\@item_info;
	}
	my @sample_info;
	push @sample_info,$s;
	for(my $i=1;$i<@{$sample_info_tmp[0]};$i++){
	    my $state=1;
	    for(my $j=0;$j<@sample_info_tmp;$j++){
		if($sample_info_tmp[$j][$i]<$state){
		    $state=$sample_info_tmp[$j][$i];
		}
	    }
	    push @sample_info,$state;
	}
	push @quality_sample,\@sample_info;
    }
    
#output qualities
    my $qf_sample=$out_dir."sample.summary";
    my $qf_file=$out_dir."file.summary";
    output_quality_data(\@header,\@quality_sample,$qf_sample);
    output_quality_data(\@header,\@quality_file,$qf_file);
    print STDERR "Merge statistics ...\n";
    my $p_sample=$out_dir."sample";
    my $p_file=$out_dir."file";
    print STDERR "Plot statistics ...\n";
    plotQual($qf_sample,$p_sample);
    plotQual($qf_file,$p_file);

    sub output_quality_data{
	my ($header_ad,$data_ad,$outfile)=@_;
	open(OUT,">$outfile")||die("Error: Unable to output to $outfile\n");
	print OUT "Sample";
	foreach my $h (@{$header_ad}){
	    print OUT "\t$h";
	}
	print OUT "\n";
	foreach my $s (@{$data_ad}){
	    for(my $i=0;$i<@{$s}-1;$i++){
		print OUT "${$s}[$i]\t";
	    }
	    print OUT "${$s}[@{$s}-1]\n";
	}
	close OUT;
    }
    sub plotQual{
        my ($file,$prefix)=@_;
	my $exec="plotQual.R";
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
          my $com="Rscript $exec $file $prefix 1>/dev/null 2>/dev/null\n";
          system($com);
          return 1;
       }
       else{
	   print STDERR ("plotQual.R cannot be found in your PATH! Omit plotting!\n");
	   return 0;
       }
    }
    1;
}

1;
