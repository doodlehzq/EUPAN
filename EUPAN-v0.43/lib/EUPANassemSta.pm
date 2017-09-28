#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-8-5
package assemStat;
sub runQuast{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_m $opt_s $opt_g);
getopts("hm:t:sg");

my $usage="\nUsage: eupan assemSta [options]  <assembly_directory> <output_directory> <QUAST_directory> <reference.fa>

eupan assemSta is used to map assembled contigs to reference and to check the statistics of assembled contigs (or scaffolds).


The script will call QUAST program, so the directory where quast.py locates is needed. 

Necessary input description:

  assembly_directory      <string>    This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, assembly results, including 
                                      files *.scafSeq and *.contig, should exist.

  output_directory        <string>    Results will be output to this directory.To avoid 
                                      overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  QUAST_directory         <string>    QUAST directory where quast.py locates.

  reference.fa            <string>    Reference sequence file (.fa or .fa.gz).

Options:
     -h                               Print this usage page.

     -t                   <int>       Threads used.
                                      Default: 1

     -m                   <int>       Minimum contig length used for assessment.
                                      Default: 500

     -g                               Check the statistics of gap-closed assemblies if -g is 
                                      enabled. In the assembly directory of each sample, 
                                      *_gc.scafSeq and *_gc.contig should exist.
                                      Default: check statistics of raw assemblies

     -s                               Check the statistics of assembled scaffolds if -s is enabled.
                                      Default: check statistics of assembled contigs

";

die $usage if @ARGV!=4;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$quast_dir,$ref)=@ARGV;


#Detect executable quast.py

$quast_dir.="/" unless($quast_dir=~/\/$/);
my $exec_quast=$quast_dir."quast.py";
die("Error01: Cannot find quast.py file in directory $quast_dir\n") unless(-e $exec_quast);

#Check existence of reference sequence file
die("Error02: Cannot find reference sequence file\n") unless(-e $ref);

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

#Read min length for assessment
my $min_length=500;
$min_length=$opt_m if(defined($opt_m));

my $suffix=".contig";
$suffix=".gcScafSeq" if(defined $opt_s && defined $opt_g);
$suffix=".scafSeq" if(!defined $opt_s && defined $opt_g);
$suffix=".gcContig" if(defined $opt_s && !defined $opt_g);


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
print STDERR "Process sample $sd ...\n"; 
#obtain *_gc.contig file within the sample directory
    my $contig_file="";
    opendir(ASS,$sd) || die("Error04: cannot open data directory: $sd\n");
    my @ass_files=readdir(ASS);
    closedir ASS;
    foreach my $f (@ass_files){
	if($f=~/$suffix$/){
	    $contig_file=$f;
	    last;
	}
    }
    if($contig_file eq ""){
	print STDERR "Warnings: cannot find assembly file(*$suffix) in $sd\n";
	next; 
    }
    $contig_file=$sd.$contig_file;

#create output directory for a sample
    my $sample_out=$out_data.$s;
    mkdir($sample_out) unless(-e $sample_out);
    my $sample_err=$err_data.$s.".err";
#generate command
    my $com;
    $com="python $exec_quast --eukaryote -t $thread_num --min-contig $min_length -o $sample_out --no-plots -R $ref $contig_file 2>&1 1>$sample_err";
    system($com);
}
my $t_sta=$out_dir."total_statistics.txt";
mergeAssemSta($out_data,$t_sta);
1;
}


sub mergeAssemSta{
    use strict;
    use warnings;
    my @input_dir=@_;
    my $outfile=pop @input_dir;
#check input file existence
    foreach my $f (@input_dir){
	die("Error01: cannot find directory $f\n") unless(-d $f);
    }
    open(OUT,">$outfile");
#read in statistics
    my %sta;
    foreach my $f (@input_dir){
	$f.="/" if $f=~/\/$/;
	opendir(DIR,$f);
	my @samples=readdir(DIR);
	closedir DIR;

	foreach my $s (@samples){
	    next if $s=~/^\.+$/;
	    my $sf=$f.$s."/report.tsv";
	    $sta{$s}=add_sta($sf) unless defined($sta{$f});
	}
    }

#output
    my $header="Assembly\t\# contigs (>= 4 bp)\tTotal length (>= 4 bp)\t\# contigs\tLargest contig\tTotal length\tReference length\tGC (%)\tReference GC (%)\tN50\tNG50\tN75\tNG75\tL50\tLG50\tL75\tLG75\t\# misassemblies\t\# misassembled contigs\tMisassembled contigs length\t\# local misassemblies\t\# unaligned contigs\tUnaligned length\tGenome fraction (%)\tDuplication ratio\t\# N's per 100 kbp\t\# mismatches per 100 kbp\t\# indels per 100 kbp\tLargest alignment\tNA50\tNGA50\tNA75\tLA50\tLGA50\tLA75";
    print OUT "Sample\t",$header,"\n";
    my @headers=split /\t/,$header;

    foreach my $s (sort keys(%sta)){
	print OUT $s;
	foreach my $h (@headers){
	print OUT "\t",$sta{$s}->{$h} if defined(${$sta{$s}}{$h});
	print OUT "\tNA" unless defined($sta{$s}->{$h});
	}
	print OUT "\n";
    }
    close OUT;
############# sub routines ##############
    sub add_sta{
	my ($file)=@_;
	print STDERR "Warning: $file doesn't exist! Not processing.\n" unless -e $file;
	return 1 unless -e $file; 
	open(FILE,$file);
	my %h;
	while(<FILE>){
	    chomp;
	    my @t=split /\t/,$_;
	    $h{$t[0]}=$t[1];
	}
	close FILE;
	return \%h;
    }
}


1;
