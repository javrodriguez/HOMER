#!/usr/bin/env perl
use warnings;


my $pthresh = 0.01;

use POSIX;

sub printCMD {
	print STDERR "\n\tUsage: cluster2bed.pl <cluster file> [options]\n";
	print STDERR "\n\tOutput: BED file to stdout\n";
	print STDERR "\n\toptions:\n";
	print STDERR "\t\t-res <#> (resolution used to create the file, default: auto detect)\n";
	print STDERR "\t\t-minp <#> (minimum % of regions in a cluster to include, default: $pthresh)\n";
	print STDERR "\t\t\t(i.e. do not output clusters containing fewer than x percent of the data\n";
	print STDERR "\t\t-name <track name> (track name for UCSC)\n";
	print STDERR "\t  i.e. cluster2bed.pl out.clusters.txt -res 1000000 -minp 0.05 > out.bed\n\n";
	exit;
}
if (@ARGV < 1) {
	printCMD();
}

my $inputFile = $ARGV[0];
my $name= "Cluster Domains from $inputFile";
my $resolution = '';

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-res') {
		$resolution = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-name') {
		$name = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-minp') {
		$pthresh = $ARGV[++$i];
	} else {
		printCMD();
	}
}

print "track name=\"$name\" itemRgb=\"On\"\n";

my %chrmax = ();
my $max = 0;
my %counts = ();
my %chr = ();
open IN, $inputFile;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$line[0] =~ /(.*?)\-(\d+)$/;
	my $chr = $1;
	my $p = $2;
	if ($p > $max) {
		$max = $p;
	}
	if (!exists($chrmax{$chr})) {
		$chrmax{$chr}=$p;
		my @a = ();
		$chr{$chr} = \@a;
	}
	if ($p > $chrmax{$chr}) {
		$chrmax{$chr}=$p;
	}
	push(@{$chr{$chr}}, $p);
	
	my $c = 0;
	if (@line > 1) {
		$c = $line[1];
		if ($c eq '') {
			$c = 0;
		}
	}
	$counts{$c}++;
	$total++;
}
close IN;

if ($resolution eq '') {
	my @diff = ();
	my $N = 0;
	foreach(keys %chr) {
		my $c = $_;
		my @pos = sort {$a <=> $b} @{$chr{$c}};
		for (my $i=1;$i<@pos;$i++) {
			push(@diff, $pos[$i]-$pos[$i-1]);
			$N++;
		}
	}
	if ($N > 0) {
		$resolution = $diff[floor(scalar(@diff)/2)];
		print STDERR "\tResulution autodetection = $resolution\n";
	} else {
		print STDERR "!!! Warning - no chromosome with more than one region - can't autodetect resolution\n";
	}
}
		

my %bad = ();
my $good = 0;
foreach(keys %counts) {
	if ($counts{$_}/$total < $pthresh) {
		$bad{$_} =1;
	} else {
		$good++;
	}
}

print STDERR "\t$good of $total clusters have more than $pthresh of data\n";

my %colors = ();
open IN, $inputFile;
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	$line[0] =~ /(.*?)\-(\d+)$/;
	my $chr = $1;
	my $p = $2;
	next if ($p == $max);
	my $p2 = $p+$resolution;
	if ($p2 > $chrmax{$chr}) {
		$p2 = $chrmax{$chr};
		next;
	}
	my $c = 0;
	if (@line > 1) {
		$c = $line[1];
		if ($c eq '') {
			$c = 0;
		}
	}
	next if (exists($bad{$c}));
	if (!exists($colors{$c})) {
		$colors{$c} = floor(rand()*255) . "," . floor(rand()*255) . "," . floor(rand()*255);
	}
	print "$chr\t$p\t$p2\t$line[0]\t0\t+\t$p\t$p2\t$colors{$c}\n";


}
close IN;
