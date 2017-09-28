#Created by Hu Zhiqiang, 2014-12-5
package geneCov;

sub cov{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_m $opt_t $opt_q);
getopts("hm:t:q:");

my $usage="\nUsage: eupanLSF geneCov [options]  <bam_directory> <output_directory> <gene_annotation>

eupanLSF geneCov is used to calculate gene coverages of each gene.

The script will call samtools and ccov.

Necessary input description:

  bam_directory           <string>    This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, mapping result, a sorted .bam
                                      file, should exist.

  output_directory        <string>    Results will be output to this directory.To avoid 
                                      overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  gene_annotation         <string>    gene annotations in a single gtf file

Options:
     -h                               Print this usage page.

     -t                   <int>       Thread number.
                                      Default:1

     -q                   <string>    LSF queue name
                                      Default: default queue    
";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$gtf)=@ARGV;

#Detect executable samtools and ccov
my $exe_ccov="ccov";
my @path=split /:/,$ENV{PATH};
my $fpflag=0;
foreach my $p (@path){
  $p.="/".$exe_ccov;
  if(-e $p && -x $p){
     $fpflag=1;
	last;
  }
}
die("ccov cannot be found in your PATH!\n
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
$thread_num=$opt_t if defined $opt_t;

#Adjust directory names and create output directory
$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

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
opendir(DATA,$data_dir) || die("Error03: cannot open data directory: $data_dir\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $sd=$data_dir.$s."/";
    next unless(-d $sd);

#obtain *.bam file within the sample directory
    my $bam_file=$sd.$s.".bam";
    unless(-e $bam_file){
	print STDERR "Warnings: cannot find bam file($bam_file) in $sd: skip this sample\n";
	next; 
    }

#create output directory for a sample
    my $sample_out=$out_data.$s;
    mkdir($sample_out) unless(-e $sample_out);
    $sample_out.="/";
    my $sta_file=$sample_out.$s.".sta";
#generate command
    my $com="$exe_ccov  $gtf $bam_file > $sta_file\n";

#generate and submit job script
#************** Might be modified for different task submission system *******************
    my $job_file=$script_out."/".$s.".lsf";   #script_file
    my $err_file=$stderr_out."/".$s.".err";   #stderr_output_file
    my $out_file=$stdout_out."/".$s.".out";   #stdout_output_file
    #create job script
    open(JOB,">$job_file")||die("Error05: Unable to create job file: $job_file\n");
    print JOB "\#BSUB -J $s","_genesta\n";              #job name
    print JOB "\#BSUB -q $opt_q\n" if defined $opt_q;  #queue name in the submission system
    print JOB "\#BSUB -o $out_file\n";                 #stdout
    print JOB "\#BSUB -e $err_file\n";                 #stderr
    print JOB "\#BSUB -n $thread_num\n";               #thread number
    print JOB "\#BSUB -R \"span[ptile=$thread_num]\"\n";

    print JOB "$com\n";                                #commands
    close JOB;
    system("bsub <$job_file");                         #submit job
#*****************************************************************************************
}
1;
}

sub merge{
    my $usage="\nUsage: eupanLSF mergeGeneCov <output_prefix> <data_directory> 
";
    die $usage if @ARGV!=2; 
    my ($prefix,$d)=@ARGV;

    my $cds_out=$prefix."cds.cov";
    my $gene_out=$prefix."gene.cov";
    my %cds;
    my %gene;
    my @genename;
    my $gf=1;

    die "invalid directory $d\n" unless -d $d;
    $d.="/" unless $d=~/\/$/;
    opendir(DIR,$d);
    my @files=readdir(DIR);
    closedir DIR;

    foreach my $f (@files){
	next if $f=~/^\.+$/;
	my ($s)=split /\./,$f;
	my $fname=$d.$f."/".$s.".sta";
	process($s,$fname) if @genename>0;
	process_and_read_name($s,$fname) if @genename==0;
    }


    my @sample=sort(keys(%gene));
    my $number=@genename;
    open(GOUT,">$gene_out");
    open(COUT,">$cds_out");
    print GOUT "Gene\t",join("\t",@sample),"\n";
    print COUT "Gene\t",join("\t",@sample),"\n";
    for(my $i=0;$i<$number;$i++){
	print GOUT $genename[$i];
	print COUT $genename[$i];
	for(my $j=0;$j<@sample;$j++){
	    print GOUT "\t",$gene{$sample[$j]}->[$i];
	    print COUT "\t",$cds{$sample[$j]}->[$i];
	}
	print GOUT "\n";
	print COUT "\n";
    }

    sub process{
	open(IN,$_[1]);
	my @ccov;
	my @gcov;
	while(<IN>){
	    next if $_=~/^\#/;
	    my @t=split /\t/,$_;
	    push @ccov,$t[7];
	    push @gcov,$t[4];
	}
	close IN;
	$gene{$_[0]}=\@gcov;
	$cds{$_[0]}=\@ccov;
    }


    sub process_and_read_name{
	open(IN,$_[1]);
	my @ccov;
	my @gcov;
	while(<IN>){
	    next if $_=~/^\#/;
	    my @t=split /\t/,$_;
	    push @genename,$t[0];
	    push @ccov,$t[7];
	    push @gcov,$t[4];
	}
	close IN;
	$gene{$_[0]}=\@gcov;
	$cds{$_[0]}=\@ccov;
    }
}
1;

