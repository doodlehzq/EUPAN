#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-10-5
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_t $opt_g $opt_s $opt_r $opt_w $opt_u $opt_c $opt_n $opt_k $opt_m);
getopts("ht:g:s:r:w:u:c:n:k:m:");

my $usage="Usage: $0 [options]  <data_directory> <output_directory> <soap_denovo2_directory>

$0 is to assembly paired FASTQ reads with SOAP DENOVO 2 with a dynamic K-mer.

ATTENTION: Executable 1)SOAPdenovo-63mer, 2)SOAPdenovo-127mer and 3)GapCloser should EXIST in soap_denovo2_directory!

Necessary input description:

  data_directory          <string>    This directory should contain one or more pair of FASTQ files
                                      with suffix of .fq or .fq.gz. The suffix can be changed with -s option.

  output_directory        <string>    The output directory.

  soap_denovo2_directory  <string>    Directory where soapdenovo2 program locates. SOAPdenovo-63mer, 
                                      SOAPdenovo-127mer and GapCloser should exist in this directory.

Options:
     -h                               Print this usage page.

     -t                   <int>       Threads used.
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

=======================================================
Any comments or bugs could be sent to Hu Zhiqiang:
    doodlehzq\@sjtu.edu.cn.
=======================================================
";

die $usage if @ARGV!=3;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$tool_dir)=@ARGV;

#check existence of output directory
#    die("Error: output directory \"$out_dir\" already exists.
#To avoid overwriting of existing files. We kindly request that the
# output directory should not exist."); if -e $out_dir;

#detect executables
$tool_dir.="/" unless $tool_dir=~/\/$/;
my $exec63=$tool_dir."SOAPdenovo-63mer";
die("Cannot find $exec63 in directory $tool_dir\n") unless -e $exec63;
my $exec127=$tool_dir."SOAPdenovo-127mer";
die("Cannot find $exec127 in directory $tool_dir\n") unless -e $exec127;
my $execgap=$tool_dir."GapCloser";
die("Cannot find $execgap in directory $tool_dir\n") unless -e $execgap;

#define thread number
my $thread_num=1;
$thread_num=$opt_t if defined($opt_t);

#define file suffix
my $suffix=".fq.gz";
$suffix=$opt_s if defined($opt_s);

#define genome size
my $gsize=380000000;
$gsize=$opt_g if defined($opt_g);

#define linear function of Kmer
my ($lin_a,$lin_b)=(0.76,20);
($lin_a,$lin_b)=split /,/,$opt_r if defined($opt_r);

#define step of Kmer change
my $step=2;
$step=$opt_w if defined($opt_w);

#define max iteration times
my $max_iter=10;
$max_iter=$opt_u if defined($opt_u);

#define SOAPdenovo contig parameters
my ($max_rd_len,$avg_ins,$reverse_seq,$asm_flags,$rd_len_cutoff,$rank,$pair_num_cutoff,$map_len)
    =(80,460,0,3,80,1,3,32);
($max_rd_len,$avg_ins,$reverse_seq,$asm_flags,$rd_len_cutoff,$rank,$pair_num_cutoff,$map_len)
    =split /,/,$opt_c if defined($opt_c);

#define the min length of contigs 
my $min_contig_len=100;
$min_contig_len=$opt_n if defined($opt_n);

#define Kmer range
my ($minKmer, $maxKmer)=(15,127); 
($minKmer, $maxKmer)=split /,/,$opt_k if defined($opt_k);

#define number of consecutive Ns
my $con_N=10;
$con_N=$opt_m if defined($opt_m);

#adjust directory and create output directory
$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);
mkdir($out_dir);

#read in files under data_directory
opendir(DATA,$data_dir) || die("Error01: can not read input data directory: $data_dir\n");
my @files=readdir(DATA);
closedir DATA;

#parse fastq pairs
my %fq_base;
foreach my $f (@files){
    next if $f=~/^\.+$/;
    print STDERR "Warnig: $f without suffix: $suffix\n" unless $f=~/$suffix$/;
    next unless $f=~/$suffix$/;
    my $fb=substr($f,0,length($f)-length($suffix)-1);
    $fq_base{$fb}=1 unless defined($fq_base{$fb});
}

#generate SOAPdenovo config
my $config_file=$out_dir."soap.config";
open(CONFIG,">$config_file") || die("Error02: can not write soap config file:$config_file\n");
print CONFIG "\#maximal read length\nmax_rd_len=$max_rd_len\n[LIB]\n\#average insert size\navg_ins=$avg_ins\n\#if sequence needs to be reversed\nreverse_seq=$reverse_seq\n\#in which part(s) the reads are used\nasm_flags=$asm_flags\n\#use only first 100 bps of each read\nrd_len_cutoff=$rd_len_cutoff\n\#in which order the reads are used while scaffolding\nrank=$rank\n\# cutoff of pair number for a reliable connection (at least 3 for short insert size)\npair_num_cutoff=$pair_num_cutoff\n\#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\nmap_len=$map_len\n\#a pair of fastq file, read 1 file should always be followed by read 2 file\n";
foreach my $b (keys(%fq_base)){
    my $forward=$data_dir.$b."1".$suffix;
    my $reverse=$data_dir.$b."2".$suffix;
    print STDERR "Warning: missed file: $forward\n" unless -e $forward;
    print STDERR "Warning: missed file: $reverse\n" unless -e $reverse;
    print STDERR "Warning: missed file: $reverse\n" unless(-e $forward && -e $reverse);
    print CONFIG "q1=$forward\nq2=$reverse\n" if(-e $forward && -e $reverse);
}
close CONFIG;

#calculate sequencing depth
my $base_number=0;
foreach my $b (keys(%fq_base)){
    my $forward=$data_dir.$b."1".$suffix;
    my $reverse=$data_dir.$b."2".$suffix;
    $base_number+=count_base_from_fq($forward);
    $base_number+=count_base_from_fq($reverse);
}
my $depth=$base_number/$gsize;
print STDERR "Total base number: $base_number\n";
print STDERR "Sequencing depth = $depth (Referense size=$gsize)\n";

#calculate init Kmer
my $Kmid=2*int(0.5*($lin_a*$depth+$lin_b))+1;
my $Klow=$Kmid-$step;
my $Khigh=$Kmid+$step;
print STDERR "Init Kmer = $Kmid\n";

#init assembly
my $Nmid=run_assembly($Kmid,$out_dir."K".$Kmid);
my $Nlow=run_assembly($Klow,$out_dir."K".$Klow);
my $Nhigh=run_assembly($Khigh,$out_dir."K".$Khigh);

#interated assembly
my $iter_times=0;
while($Nmid<$Nlow || $Nmid<$Nhigh){
    $iter_times++;
    if($Nmid<$Nhigh){        #Kmer goes higher
	if($iter_times>$max_iter){        #check iteration times
	    print STDERR "Stop best Kmer selection: out of max iteration times($max_iter)\n";
	    ($Kmid,$Khigh)=($Khigh,$Kmid);
	    ($Nmid,$Nhigh)=($Nhigh,$Nmid);
	    last;
	}
	if($Khigh+$step>$maxKmer){        #check Kmer region
	    print STDERR "Stop best Kmer selection: ouside the UP boundary of Kmer($maxKmer)\n";
	    ($Kmid,$Khigh)=($Khigh,$Kmid);
	    ($Nmid,$Nhigh)=($Nhigh,$Nmid);
	    last;
	}
	$Klow=$Kmid;
	$Nlow=$Nmid;
	$Kmid=$Khigh;
	$Nmid=$Nhigh;
	$Khigh+=$step;
	$Nhigh=run_assembly($Khigh,$out_dir."K".$Khigh);
    }
    elsif($Nmid<$Nlow){        #Kmer goes lower
	if($iter_times>$max_iter){        #check iteration times
	    print STDERR "Stop best Kmer selection: out of max iteration times($max_iter)\n";
	    ($Kmid,$Klow)=($Klow,$Kmid);
	    ($Nmid,$Nlow)=($Nlow,$Nmid);
	    last;
	}
	if($Klow-$step<$minKmer){        #check Kmer region
	    print STDERR "Stop best Kmer selection: ouside the DOWN boundary of Kmer($minKmer)\n";
	    ($Kmid,$Klow)=($Klow,$Kmid);
	    ($Nmid,$Nlow)=($Nlow,$Nmid);
	    last;
	}
	$Klow=$Kmid;
	$Nlow=$Nmid;
	$Kmid=$Khigh;
	$Nmid=$Nhigh;
	$Khigh+=$step;
	$Nhigh=run_assembly($Khigh,$out_dir."K".$Khigh);
    }
}

#adjust outputs
system("mv $out_dir"."K$Kmid/* $out_dir");                                                                                                                                     
print STDERR "Finished!\nIteration_time=$iter_times\nKmer(selected)=$Kmid\nN50(K=$Kmid)=$Nmid\n";
exit;



######################## sub routines ############################

sub count_base_from_fq{      #usage: count_base_from_fq(file_name)                                                                                                                  
    my ($fq)=@_;
    my $bn=0;
    if($fq=~/\.gz$/){
        open(FQ,"zcat $fq |") ||die("Error07: cannot open gzipped file: $fq\n");
    }
    else{
        open(FQ,"$fq") ||die("Error08: cannot open file: $fq\n");
    }
    my $i=0;
    while(<FQ>){
        $i++;
        chomp;                                                                                                                                                                    
        $bn+=length($_) if($i%4==2);                                                                                                                                             
    }
    close FQ;
    return $bn;
}


sub run_assembly{          #usage n50=run_assembly(k,outdir)
    my ($kmer,$out)=@_;
    print STDERR "Assembly with kmer=$kmer\n";
    mkdir($out) unless(-e $out);

#run soap_denovo                                                                                                                                                                  
    my $com;
    if($kmer<=63){
        $com=$exec63;
    }
    else{
	$com=$exec127;
    }
    $com.=" all -s $config_file -o $out/K$kmer -K $kmer -R -F -p $thread_num >$out/soap.log 2>&1";
    system($com);

#run gapcloser
    $com=$execgap;
    $com.=" -b $config_file -a $out/K$kmer.scafSeq -o $out/K\Q$kmer\E.gcScafSeq -t $thread_num >$out/gapcloser.log 2>&1";
    system($com);

#break gap closed scaffolds to contigs
    break_scaffolds("$out/K\Q$kmer\E.gcScafSeq","$out/K\Q$kmer\E.gcContig",$con_N);

#calculate N50

    my $n50=calculate_N50("$out/K\Q$kmer\E.gcContig");
    print STDERR "N50(K=$kmer)=$n50\n";
    return $n50;
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

sub calculate_N50{  #usage: n50=calculate_N50(contig_file)
    open(CONTIG,$_[0])||die("Error11: cannot open file:$_[0]\n");
    my @len;
    my $seq="";
    my $cl=0;
    my $init=1;
    while(<CONTIG>){
        if(/^>/){
            if(!$init){
                push @len,$cl if($cl>=$min_contig_len);
                $seq="";
                $cl=0;
            }
            else{
                $init=0;
            }
        }
        else{
            chomp;
            $cl+=length($_);
        }
    }
    push @len,$cl if($cl>=$min_contig_len);
    close CONTIG;
    @len=sort{$a<=>$b}(@len);
    my $total=0;
    foreach my $v (@len){
        $total+=$v;
    }
    $total/=2;
    my $cal=0;
    for(my $i=0;$i<@len;$i++){
        $cal+=$len[$i];
        return $len[$i] if($cal>=$total);
    }
}
