#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-9-5
package assembly;

sub assemble{
    my $usage="
Usage: eupan assemble [commands] ...

Commands:
\tsoapdenovo    Aseembly with SOAPdenovo 2.
\tlinearK       Assembly with an iterative use of SOAPdenovo 2 (Recommended).
";
    die $usage if @ARGV<1;
    my $com=shift @ARGV;
    if($com eq "soapdenovo"){
	soap(@ARGV);
    }
    elsif($com eq "linearK"){
	linearK(@ARGV);
    }
    else{
	print STDERR "Unknown command: $com\n";
	die($usage);
    }
}

sub soap{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_s $opt_n $opt_m $opt_k $opt_s $opt_c $opt_g);
getopts("hs:n:m:t:k:s:c:g");

my $usage="\nUsage: eupan assemble soapdenovo [options] <fastq_data_directory> <output_directory> <soapdenovo_directory>

eupan assemble soapdenovo is used to assemble high-quality reads on large scale.

Necessary input description:

  fastq_data_directory    <string>    This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, there should be several 
                                      sequencing files ended by .fastq or .fastq.gz.

  output_directory        <string>    Alignment results will be output to this directory.
                                      To avoid overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  sopadenovo_directory    <string>    directory where soapdenovo2 executable files exists   

Options:
     -h                              Print this usage page.

     -t                   <int>      Threads used.
                                     Default: 1

     -s                   <string>    Suffix of files within data_directory.
                                      Default: .fq.gz 

     -k                   <int>      Kmer.
                                     Default: 35

     -c                   <string>    Parameters of soapdenovo2 config file. 8 parameters ligated by comma
                                        1)maximal read length
                                        2)average insert size
                                        3)if sequence needs to be reversed
                                        4)in which part(s) the reads are used
                                        5)use only first N bps of each read
                                        6)in which order the reads are used while scaffolding
                                        7)cutoff of pair number for a reliable connection (at least 3 for 
                                          short insert size)
                                        8)minimum aligned length to contigs for a reliable read location 
                                          (at least 32 for short insert size)
                                      Default: 80,460,0,3,80,1,3,32

     -g                               enable gapcloser 
";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$tool_dir)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
}

$tool_dir.="/" unless $tool_dir=~/\/$/;
my $exec63=$tool_dir."SOAPdenovo-63mer";
die("Cannot find $exec63 in directory $tool_dir\n") unless -e $exec63;
my $exec127=$tool_dir."SOAPdenovo-127mer";
die("Cannot find $exec127 in directory $tool_dir\n") unless -e $exec127;
my $execgap=$tool_dir."GapCloser";
if(defined $opt_g){
die("Cannot find $execgap in directory $tool_dir\n") unless -e $execgap;
}
#read threads
my $thread_num=1;
if(defined($opt_t)){
    $thread_num=$opt_t;
}
#define kmer
my $kmer=35;
$kmer=$opt_k if defined $opt_k;

#define file suffix
my $suffix=".fq.gz";
$suffix=$opt_s if defined($opt_s);

my ($max_rd_len,$avg_ins,$reverse_seq,$asm_flags,$rd_len_cutoff,$rank,$pair_num_cutoff,$map_len)
    =(80,460,0,3,80,1,3,32);
($max_rd_len,$avg_ins,$reverse_seq,$asm_flags,$rd_len_cutoff,$rank,$pair_num_cutoff,$map_len)
    =split /,/,$opt_c if defined($opt_c);

#define the min length of contigs 
my $min_contig_len=500;
$min_contig_len=$opt_n if defined($opt_n);


#Adjust directory names and create output directory

$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);

#read samples
opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $s_dir=$data_dir.$s."/";
    print STDERR "Process sample $s_dir ...\n";
    my $o_dir=$out_data.$s."/";
    mkdir($o_dir);
    opendir(DATA,$s_dir) || die("Error: can not read input data directory: $s_dir\n");
    my @files=readdir(DATA);
    closedir DATA;

#parse fastq pairs
    my %fq_base;
    foreach my $f (@files){
	next if $f=~/^\.+$/;
	next if $f=~/^single/;
	print STDERR "Warnig: $f without suffix: $suffix\n" unless $f=~/$suffix$/;
	next unless $f=~/$suffix$/;
	my $fb=substr($f,0,length($f)-length($suffix)-1);
	$fq_base{$fb}=1 unless defined($fq_base{$fb});
    }

#generate SOAPdenovo config
    my $config_file=$o_dir."soap.config";
    open(CONFIG,">$config_file") || die("Error02: can not write soap config file:$config_file\n");
    print CONFIG "\#maximal read length\nmax_rd_len=$max_rd_len\n[LIB]\n\#average insert size\navg_ins=$avg_ins\n\#if sequence needs to be reversed\nreverse_seq=$reverse_seq\n\#in which part(s) the reads are used\nasm_flags=$asm_flags\n\#use only first 100 bps of each read\nrd_len_cutoff=$rd_len_cutoff\n\#in which order the reads are used while scaffolding\nrank=$rank\n\# cutoff of pair number for a reliable connection (at least 3 for short insert size)\npair_num_cutoff=$pair_num_cutoff\n\#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\nmap_len=$map_len\n\#a pair of fastq file, read 1 file should always be followed by read 2 file\n";
    foreach my $b (keys(%fq_base)){
	my $forward=$s_dir.$b."1".$suffix;
	my $reverse=$s_dir.$b."2".$suffix;
	print STDERR "Warning: missed file: $forward\n" unless -e $forward;
	print STDERR "Warning: missed file: $reverse\n" unless -e $reverse;
	print STDERR "Warning: missed file: $reverse\n" unless(-e $forward && -e $reverse);
	print CONFIG "q1=$forward\nq2=$reverse\n" if(-e $forward && -e $reverse);
    }
    close CONFIG;
    my $com;
    if($kmer<=63){
        $com=$exec63;
    }
    else{
	$com=$exec127;
    }

    $com.=" all -s $config_file -o $o_dir/K$kmer -K $kmer -R -F -p $thread_num >$o_dir/soap.log 2>&1\n";
    system($com);
#run gapcloser
    if(defined $opt_g){
	$com=$execgap;
	$com.=" -b $config_file -a $o_dir/K$kmer.scafSeq -o $o_dir/K\Q$kmer\E.gcScafSeq -t $thread_num >$o_dir/gapcloser.log 2>&1";
	system($com);
	break_scaffolds("$o_dir/K\Q$kmer\E.gcScafSeq","$o_dir/K\Q$kmer\E.gcContig",10);
    }
}

sub break_scaffolds{
    my ($in,$out,$sepN)=@_;
    my $sep="";
    for(my $i=1;$i<$sepN;$i++){
	$sep.="N";
    }
    open(SCAF,$in)||die("Error09: cannot read file:$in\n");
    open(CONTIG,">$out")||die("Error10: cannot write file:$out\n");
    my $init=1;
    my ($header,$seq);
    while(my $line=<SCAF>){
	chomp $line;
	if($line=~/^>/){
	    if(!$init){
		split_and_print_seq($seq,$sep,$header);
		($header)=split /\s+/,$line;
		$seq="";
	    }
	    else{
		($header)=split /\s+/,$line;
		$seq="";
		$init=0;
	    }
	}
	else{
	    $seq.=$line;
	}
    }
    split_and_print_seq($seq,$sep,$header);
    close SCAF;
    close CONTIG;
}



sub split_and_print_seq{
    my ($seq,$sep,$header)=@_;
    my $i=0;
    my @tmp=split /\Q$sep\EN+/,$seq;
    foreach my $s (@tmp){
	next if length($s)==0;
	$i++;
	print CONTIG $header,"_",$i,"\n";
	for(my $j=0;$j<int(length($s)/100);$j++){
	    print CONTIG substr($s,100*$j,100),"\n";
	}
	print CONTIG substr($s,0-length($s)%100),"\n" if length($s)%100>0;
    }
}

1;
}

sub linearK{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_s $opt_g $opt_r $opt_w $opt_u $opt_c $opt_n $opt_k $opt_m);
getopts("ht:g:s:r:w:u:c:n:k:m:");

my $usage="\nUsage: eupan assemble linearK [options] <fastq_data_directory> <output_directory> <soapdenovo_directory>

eupan assemble linearK is used to assemble high-quality reads on large scale.

Necessary input description:

  fastq_data_directory    <string>    This directory should contain many sub-directories
                                      named by sample names, such as CX101, B152,etc.
                                      In each sub-directory, there should be several 
                                      sequencing files ended by .fastq or .fastq.gz.

  output_directory        <string>    Alignment results will be output to this directory.
                                      To avoid overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

  sopadenovo_directory    <string>    directory where soapdenovo2 executable files exists   

Options:
     -h                              Print this usage page.

     -t                   <int>      Threads used.
                                     Default: 1

     -g                   <int>       Genome size. Used to infer sequencing depth. 
                                      Default: 380000000 (460M)
     
     -s                   <string>    Suffix of files within data_directory.
                                      Default: .fq.gz 

     -r                   <string>    Parameters of linear function: Kmer=2*int(0.5*(a*Depth+b))+1. 
                                      The parameter should be input as \"a,b\".
                                      Default: 0.76,20

     -w                   <int>       Step-length of Kmer change.
                                      Default: 2

     -u                   <int>       Upper limmited times of Kmer change. This parameter is set to reduce
                                      redundancy computation.
                                      Default: 10

     -c                   <string>    Parameters of soapdenovo2 config file. 8 parameters ligated by comma
                                        1)maximal read length
                                        2)average insert size
                                        3)if sequence needs to be reversed
                                        4)in which part(s) the reads are used
                                        5)use only first N bps of each read
                                        6)in which order the reads are used while scaffolding
                                        7)cutoff of pair number for a reliable connection (at least 3 for 
                                          short insert size)
                                        8)minimum aligned length to contigs for a reliable read location 
                                          (at least 32 for short insert size)
                                      Default: 80,460,0,3,80,1,3,32

     -n                   <int>       The minimum length of contigs. Contigs shorter than this length will
                                      NOT be used when calculating N50.
                                      Default: 100

     -k                   <string>    Available Kmer range. Give comma-seperated lower bound and upper bound.
                                      Default: 15,127

     -m                   <int>       The number of consecutive Ns to be broken down to contigs.This is used 
                                      in the process break gapclosed scaffolds to contigs.
                                      Default: 10.
";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$soapdenovo)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist
.");
}
#Check executable linearK
my $exec="linearK";
my @path=split /:/,$ENV{PATH};
my $fpflag=0;
foreach my $p (@path){
  $p.="/".$exec;
  if(-e $p && -x $p){
     $fpflag=1;
	last;
  }
}
die("Executable linearK cannot be found in your PATH!\n
") unless($fpflag);


#read threads
my $thread_num=1;
if(defined($opt_t)){
    $thread_num=$opt_t;
}

#Adjust directory names and create output directory

$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);

#read samples
opendir(DATA,$data_dir) || die("Error: can not open input data directory!\n");
my @sample=readdir(DATA);
closedir DATA;

#process each sample
foreach my $s (@sample){
    next if $s=~/^\./;
    my $com="linearK -t $thread_num";
    $com.=" -g $opt_g" if defined $opt_g;
    $com.=" -s $opt_s" if defined $opt_s;
    $com.=" -r $opt_r" if defined $opt_r;
    $com.=" -w $opt_w" if defined $opt_w;
    $com.=" -u $opt_u" if defined $opt_u;
    $com.=" -c $opt_c" if defined $opt_c;
    $com.=" -n $opt_n" if defined $opt_n;
    $com.=" -k $opt_k" if defined $opt_k;
    $com.=" -m $opt_m" if defined $opt_m;
    my $indir=$data_dir.$s;
    print STDERR "Process sample $indir ...\n";
    my $outdir=$out_data.$s;
    $com.=" $indir $outdir $soapdenovo";
    system($com);                           #submit job
#*****************************************************************************************
}
1;
}

1;

