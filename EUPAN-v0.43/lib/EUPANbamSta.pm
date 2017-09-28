#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-8-5
package bamStat;
sub bamsta{
    my $usage="
eupan bamSta [commands] ...

Commands:
\tbasic               calculate basic statistics
\tcov                 calculate genome coverage
";
    die $usage if @ARGV<1;
    my $com=shift @ARGV;
    if($com eq "basic"){
	basic(@ARGV);
    }
    elsif($com eq "cov"){
	cov(@ARGV);
    }
    elsif($com eq "mergeBasicSta"){
	mergeBasicSta(@ARGV);
    }
    elsif($com eq "mergeCovSta"){
	mergeCovSta(@ARGV);
    }
    else{
	print STDERR "Unknown command: $com\n";
	die($usage);
    }
}

sub basic{
use strict;
use warnings;
my $usage="\nUsage: eupan bamSta basic [options]  <bam_directory> <output_directory> <bamUtil_directory>

eupan bamSta basic is used to check the basic statistics of mapping.

The script will call bam_stats (in BamUtil), so the directory where bamUtil locates is needed. 

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

  bamUtil_directory       <string>    bamUtil directory where bin/bam locates.


";

die $usage if @ARGV<3;
my ($data_dir,$out_dir,$bamutil_dir)=@ARGV;


#Detect executable bam_stats

$bamutil_dir.="/" unless($bamutil_dir=~/\/$/);
my $exec_bs=$bamutil_dir."bin/bam";
die("Error01: Cannot find bam_stats file in directory bin/ under $bamutil_dir\n") unless(-e $exec_bs);

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
    my $sample_out=$out_data.$s;
    mkdir($sample_out) unless(-e $sample_out);
    $sample_out.="/$s.sta";
#generate command
    my $com;
    $com="$exec_bs stats --basic --in $bam_file 2>$sample_out";
    system($com);
}
my $summary=$out_dir."summary.txt";
mergeBasicSta($out_data,$summary);
1;
}

sub cov{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_m $opt_t $opt_q);
getopts("hm:t:q:");

my $usage="\nUsage: eupan bamSta cov [options]  <bam_directory> <output_directory> <qualimap_directory>

eupan bamSta cov is used to check the coverages of the genome.

The script will call qualimap, so the directory where qualimap locates is needed. 

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

  qualimap_directory      <string>    qualimap directory where executable qualimap locates.

Options:
     -h                               Print this usage page.

     -m                   <string>    Maximum memory size for java use
                                      Default: 12G

     -t                   <int>       Thread number.
                                      Default:4

     -q                   <string>    LSF queue name
                                      Default: default queue    

";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$qualimap_dir)=@ARGV;


#Detect executable bam_stats

$qualimap_dir.="/" unless($qualimap_dir=~/\/$/);
my $exec_bs=$qualimap_dir."qualimap";
die("Error01: Cannot find qualimap file in directory bin/ under $qualimap_dir\n") unless(-e $exec_bs);

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
}

#Read threads
my $thread_num=4;
$thread_num=$opt_t if defined $opt_t;

#Max mem
my $max_mem="12G";
$max_mem=$opt_m if defined $opt_m;

#Adjust directory names and create output directory
$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

#Create output directory and sub-directories
mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);
my $err_data=$out_dir."err/";
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
#obtain *.bam file within the sample directory
    my $bam_file=$sd.$s.".bam";
    unless(-e $bam_file){
	print STDERR "Warnings: cannot find bam file($bam_file) in $sd: skip this sample\n";
	next; 
    }

#create output directory for a sample
    my $sample_out=$out_data.$s;
    mkdir($sample_out) unless(-e $sample_out);
    my $sample_err=$err_data.$s.".err";

#generate command
    my $com;
    $com="$exec_bs bamqc --bam $bam_file --outdir $sample_out --java-mem-size=$max_mem 1>$sample_err\n";
    $com.="rm $sample_out"."/coverage.txt\n";
    system($com);
}
my $summary=$out_dir."summary.txt";
mergeCovSta($out_data,$summary);
1;
}

sub mergeBasicSta{
my ($data_dir,$outfile)=@_;

my $title="Sample\tTotalReads(e6)\tMappedReads(e6)\tPairedReads(e6)\tProperPair(e6)\tDuplicateReads(e6)\tQCFailureReads(e6)\tMappingRate(%)\tPairedReads(%)\tProperPair(%)\tDupRate(%)\tQCFailRate(%)\tTotalBases(e6)\tBasesInMappedReads(e6)\n";
my %data;
$data_dir.="/" unless $data_dir=~/\/$/;
opendir(DIR,$data_dir);
my @sample=readdir(DIR);
closedir DIR;

foreach my $s (@sample){
    next if $s=~/^\.+$/;
    my $sd=$data_dir.$s."/";
    unless(-d $sd){
	print STDERR "$sd isn't a directory! Omit!\n";
	next;
    }
    opendir(DIR,$sd);
    my @file=readdir(DIR);
    closedir DIR;
    foreach my $ff (@file){
	next if $ff=~/^\.+$/;
	my $f=$sd.$ff;
	if($ff=~/^(.+)\.sta$/){
	    my $name=$1;
	    my @val=();
	    open(IN,$f);
	    while(<IN>){
		chomp;
		if(/^Number/){
		    next;
		}
		elsif(/\w+/){
		    my @t=split /\s+/,$_;
		    if(@t==2){
			push @val, $t[1];
		    }
		}
	    }
	    close IN;
	    $data{$name}=join("\t",@val);
	}
	else{
	    print STDERR "$f isn't ended with .sta\nOmit!\nPlease check the input directory structure!
";
	    next;
	}
    }
}
open(OUT,">$outfile");
print OUT $title;
foreach my $s (sort keys(%data)){
    print OUT $s,"\t",$data{$s},"\n";
}
close OUT;
1;
}

sub mergeCovSta{
    use strict;
    use warnings;
    my ($data_dir,$outfile)=@_;
    $data_dir.="/" unless $data_dir=~/\/$/;
# There is a 85.71% of reference with a coverageData >= 2X
    my $title="Sample\tCoverage(%)\n";
    my %data;
    opendir(DIR,$data_dir);
    my @sample=readdir(DIR);
    closedir DIR;

    foreach my $s (@sample){
	next if $s=~/^\.+$/;
	my $sd=$data_dir.$s."/";
	unless(-d $sd){
	    print STDERR "$sd isn't a directory! Omit!\n";
	    next;
	}
	my $f=$sd."genome_results.txt";
	unless(-e $f){
	    print STDERR "File $f doesn't exist!  Omit sample $s\n";
	    next;
	}
	open(IN,$f);
	while(<IN>){
	    chomp;
	    if(/^\s+There is a ([0-9\.]+)% of reference with a coverageData >= 2X/){
		$data{$s}=$1;
	    }
	}
	close IN;
    }
    open(OUT,">$outfile");
    print OUT $title;
    foreach my $s (sort keys(%data)){
	print OUT $s,"\t",$data{$s},"\n";
    }
    close OUT;
    1;
}
1;


