#!/usr/bin/perl
#Created by Hu Zhiqiang, 2014-7-2
package qualSta;
sub checkQual{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_f $opt_t $opt_v);
    getopts("hf:t:v:");
    my $usage="\nUsage: eupan qualSta [options] <data_directory> <output_directory>

qualSta is used to check qualities of .fastq/.fastq.gz files on a large scale.

The script will call fastqc program, so please make sure fastqc is in your
PATH, or you need to use -f option to tell the program where fastqc locates.

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
    
     -v            <string>      Sets to: PE, if all files are PE data.
				
";

    die $usage if @ARGV!=2;
    die $usage if defined($opt_h);
    my ($data_dir,$out_dir)=@ARGV;

#Check existence of output directory
    if(-e $out_dir){
	die("Error: output directory \"$out_dir\" already exists.
To avoid overwriting of existing files. We kindly request that the
 output directory should not exist.
");
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
    print STDERR "Processing samples:\n";

#process each sample
    foreach my $s (@sample){
	next if $s=~/^\./;
	my $sd=$data_dir.$s;
	if(-d $sd){
	    #read sample directories
	    opendir(RUN,$sd)|| die("Error: can not open directory: $sd\n");
	    my @files=readdir(RUN);
	    close RUN;
	    print STDERR "\tProcess sample $sd ...\n";
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
	    system($command); 
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
    use Getopt::Std;
    use vars qw($opt_h $opt_f $opt_t);
    getopts("hf:t:");
    my $out_dir=$ARGV[1];
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
	    if(defined($opt_v)){
		if($opt_v eq "PE"){
		next if $f=~/^single/;
		}
	    }
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
    print STDERR "Merge statistics ...\n";
    my $qf_sample=$out_dir."sample.summary";
    my $qf_file=$out_dir."file.summary";
    output_quality_data(\@header,\@quality_sample,$qf_sample);
    output_quality_data(\@header,\@quality_file,$qf_file);
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

