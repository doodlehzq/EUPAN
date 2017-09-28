use strict;
use warnings;
package fastaSta;

sub sta{
    my $usage="eupan fastaSta <fasta>
";
    die $usage if @ARGV!=1;

    open(IN,$ARGV[0]);
    my $i=0;
    my $j=0;
    while(<IN>){
	chomp;
	if(/^>(.+)$/){
	    $j++;
	}
	else{
	    $i+=length($_);
	}
    }

    print STDOUT "Total sequence number: ",$i,"\n";
    print STDOUT "Total base number: ",$j,"\n";
}
1;
