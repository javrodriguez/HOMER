#!/usr/bin/env perl
use warnings;
use lib "/gpfs/data01/cbenner/software/homer/.//bin";
my $homeDir = "/gpfs/data01/cbenner/software/homer/./";

# Copyright 2009-2016 Christopher Benner <cbenner@salk.edu>
# 
# This file is part of HOMER
#
# HOMER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HOMER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


use POSIX;
use HomerConfig;
use Statistics;


my $reduceThresh = 0.6;
my $percentSimilar = 0.20;
my $matchThresh = "T10";
my $knownPvalueThresh = 0.01;
my $motifInfoFileName = "motifFindingParameters.txt";
my $deleteFlag = 1;
my %toDelete = ();
my $overlapThreshold = 0.85;
my $method = 'sample';

my $config = HomerConfig::loadConfigFile();

if (@ARGV < 2) {
	printCMD();
}

sub printCMD {
	print STDERR "\n\tThis Program will export potential background regions based on repeats\n";
	print STDERR "\n\tUsage: selectRepeatBg.pl <pos file> <genome> [additional options]\n";
	print STDERR "\n\tOutput: sends position/peak file to stdout\n";

	print STDERR "\n\tPossible Genomes:\n";
	foreach(keys %{$config->{'GENOMES'}}) {
		print STDERR "\t\t$_\t$config->{'GENOMES'}->{$_}->{'org'}\n";
	}
	print STDERR "\n\tBasic options:\n";
	#print STDERR "\t\t-mask (mask repeats/lower case sequence, can also add 'r' to genome, i.e. mm9r)\n";
	print STDERR "\t\t-f <#> (Number of bg regions per peak to select)\n";
	print STDERR "\n";
	exit;
}


my $posFile = $ARGV[0];
my $genome = $ARGV[1];
my $N = 50000;
my $factor = 5;
my $topRepeatN = 20;
if ($genome =~ s/r$//) {
	$cmd->{'mask'} = 1;
}
$genome = $genome;
print STDERR "\tPosition file = $posFile\n";
print STDERR "\tGenome = $genome\n";
	
for (my $i=2;$i<@ARGV;$i++) { 
	if ($ARGV[$i] eq '-N') {
		$N = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-f') {
		$factor = $ARGV[++$i];
	} else {
		print STDERR "!!! $ARGV[$i] not recognized!\n";
		printCMD();
	}
}	

my $genomeDir = "";
my $preparsedDirFromConfig = "";
my $customGenome = "";
my $organism = 'null';
if (!exists($config->{'GENOMES'}->{$genome})) {
	$customGenome = $genome;
	($genome,$genomeDir,$preparsedDirFromConfig) = HomerConfig::parseCustomGenome($genome);
} else {
	$genomeDir = $config->{'GENOMES'}->{$genome}->{'directory'} . "/";
	$organism =  $config->{'GENOMES'}->{$genome}->{'org'};
	$preparsedDirFromConfig = $genomeDir . "preparsed/";
}

open IN, $posFile or die "!!! Could not open peak/position file $cmd->{'posfile'} !!!\n";
close IN;


my $tmpID = rand();
my $tmpFile = $tmpID . ".tmp";
my $tmpFile2 = $tmpID . ".2.tmp";

my $mflag = "";
$mflag = " -mask " if ($cmd->{'mask'});

print STDERR "\tReading peak file...\n";
my %hits = ();
my %peaks = ();
my %peakInfo = ();
my $numPeaks = 0;
my $avgSize = 0;

open IN, $posFile or die "!!! Could not open peak file ($posFile)\n";
while (<IN>) {
	chomp;
	s/\r//g;
	my @line = split /\t/;
	next if (/^\s*#/);
	next unless ($line[2] =~ /^\d+$/);
	next unless ($line[3] =~ /^\d+$/);
	my $id = $line[0];
	my $chr = $line[1];
	my $start = $line[2];
	my $end = $line[3];
	my $strand = $line[4];
	if ($strand eq '1') {
		$strand = '-';
	} elsif ($strand eq '0') {
		$strand = '+';
	}

	$numPeaks++;
	$avgSize += $end-$start;
	my $p = {id=>$id,c=>$chr,s=>$start,e=>$end,d=>$strand,r=>''};

	if (!exists($peaks{$chr})) {
		my @a = ();
		$peaks{$chr} = \@a;
	}
	push(@{$peaks{$chr}}, $p);
	$peakInfo{$id} = $p;

}
close IN;
if ($numPeaks < 1) {
	print STDERR "!!! No peaks found in peak file\n";
	exit;
}
$avgSize = floor($avgSize/$numPeaks);
print STDERR "\t\t$numPeaks peaks found, avg size = $avgSize\n";



my $repeatFile = $genomeDir . "/$genome.repeats";
print STDERR "\tReading repeat file ($repeatFile)...\n";

my %repeats = ();
my %repInfo = ();
open IN, $repeatFile or die "!!! Could not open repeat info file: \"$repeatFile\"\n";
my $totalRepeats = 0;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^\s*#/);
	my @line = split /\t/;
	$totalRepeats++;
	my $div = $line[6]/1000.0;
	my $chr = $line[1];
	my $start = $line[2]+1;
	my $end = $line[3];
	my $strand = $line[4];
	if ($strand eq '1') {
		$strand = '-';
	} elsif ($strand eq '0') {
		$strand = '+';
	}

	my $ogStart = $line[7];
	my $ogEnd = $line[8];
	my $ogLength  = 1;
	if ($ogStart =~ /^[\d\.\-\e\+]+$/ && $ogEnd =~ /^[\d\.\-\e\+]+$/) {
		$ogLength  = $ogEnd-$ogStart;
	}
	$ogLength = 1 if ($ogLength < 1);
	my $length = $end - $start;
	my $fracLength = ($length)/$ogLength;

	my $name = $line[0];
	$name =~ s/\-HOMER.*$//;
	my $sname = $name . "-$chr:$start-$end";
	my $rand = rand();
	
	my $r = {s=>$start,e=>$end,ogs=>$ogStart,oge=>$ogEnd,c=>$chr,d=>$strand,div=>$div,
			len=>$length,oglen=>$ogLength,frac=>$fracLength,n=>$name,mask=>0,sel=>0,rand=>$rand};

	if (!exists($repInfo{$name})) {
		my @a = ();
		$repInfo{$name}=\@a;
	}
	push(@{$repInfo{$name}},$r);

	if (!exists($repeats{$chr})) {
		my @a = ();
		$repeats{$chr} = \@a;
	}
	push(@{$repeats{$chr}}, $r);
}
close IN;
print STDERR "\t\t$totalRepeats total repeats found\n";

print STDERR "\tAssigning peaks to repeat regions...\n";
foreach(keys %peaks) {
	my $chr = $_;
	next if (!exists($repeats{$chr}));
	my @peaks = sort {$a->{'s'} <=> $b->{'s'}} @{$peaks{$chr}};
	my @repeats = sort {$a->{'s'} <=> $b->{'s'}} @{$repeats{$chr}};

	#populate non-repeat genome entries
	for (my $i=0;$i<@repeats-1;$i++) {
		my $span = $repeats[$i+1]->{'s'} - $repeats[$i]->{'e'};
		my $start = $repeats[$i]->{'e'};
		my $name= 'none';
		next if ($span < $avgSize);
		for (my $j=0;$j<$span;$j+=$avgSize) {
			my $s = $start+$j;
			my $e = $s+$avgSize;
			my $rand = rand();
			my $r = {s=>$s,e=>$e,ogs=>$s,oge=>$e,c=>$chr,d=>'+',div=>0,
				len=>$avgSize,oglen=>$avgSize,frac=>1,n=>$name,mask=>0,sel=>0,score=>0,rand=>$rand};
			#my $id = $name . '-' . $chr . "-" . $i . "-" . $j;
			if (!exists($repInfo{$name})) {
				my @a = ();
				$repInfo{$name}=\@a;
			}
			push(@{$repInfo{$name}},$r);
		}
	}



	my $rindex = 0;
	for (my $pindex=0;$pindex<@peaks;$pindex++) {
		
		my $pstart = $peaks[$pindex]->{'s'};
		my $pend = $peaks[$pindex]->{'e'};
		my $pstrand = $peaks[$pindex]->{'d'};

		my $bestCov = -1;
		my $best = {cov=>0, n=>'none',coords=>'',div=>0};

		for (my $j=$rindex;$j<@repeats;$j++) {
			#look for overlaping repeats...
			my $rstart = $repeats[$j]->{'s'};
			my $rend = $repeats[$j]->{'e'};
			my $rstrand = $repeats[$j]->{'d'};
			if ($rend < $pstart) {
				$rindex++ if ($j== $rindex);
				next;
			}
			last if ($rstart > $pend);
			if (($rstart >= $pstart && $rstart <= $pend) ||
					($rend >= $pstart && $rend <= $pend) || 
					($rstart < $pstart && $rend > $pend)) {
				#Overlap!
				#get repeat region that is effected relative to mid-point
				my $ostart = $rstart;
				if ($rstart < $pstart) {
					$ostart = $pstart;
				}
				my $oend = $pend;
				if ($rend < $pend) {
					$oend = $rend;
				}
				my $offsetStart =0;
				my $offsetEnd =0;
				my $ogstart = $repeats[$j]->{'ogs'};
				my $ogend = $repeats[$j]->{'oge'};
				if ($rstrand eq '+') {
					$offsetStart = $ostart-$ogstart;
					$offsetEnd = $oend-$ogstart;
				} else {
					$offsetStart = $ogend-$oend;
					$offsetEnd = $ogend-$ostart;
				}
				$repeats[$j]->{'mask'} = 1;
				my $name = $repeats[$j]->{'n'};
				my @coords = ($offsetStart, $offsetEnd);
				my $cov=$offsetEnd-$offsetStart;
				my $div = $repeats[$j]->{'div'};
				my $x = {cov=>$cov, n=>$name,coords=>\@coords,div=>$div};
				if ($cov > $bestCov) {
					$bestCov = $cov;
					$best = $x;
				}
			}
		}
		if ($bestCov > 0) {
			my $name = $best->{'n'};
			if (!exists($hits{$name})) {
				my @a = ();
				my @b = ();
				my $x = {n=>0,c=>\@a,d=>\@b};
				$hits{$name} = $x;
			}
			$hits{$name}->{'n'}++;
			push(@{$hits{$name}->{'c'}},$best->{'coords'});
			push(@{$hits{$name}->{'d'}},$best->{'div'});
		}
		$peaks[$pindex]->{'r'} = $best;
	}
}


print STDERR "\tShowing top $topRepeatN repeat families in peaks:\n";
my @a = sort {$hits{$b}->{'n'} <=> $hits{$a}->{'n'}} keys %hits;
my $c = 0;
foreach(@a) {
	$c++;
	last if ($c > $topRepeatN);
	my $name = $_;
	my $n = $hits{$name}->{'n'};
	my $frac = sprintf("%.2f%%", $n/$numPeaks*100);
	my $avgpos = 0;
	my $avgcov = 0;
	my $avgdiv = 0;
	foreach(@{$hits{$name}->{'c'}}) {
		my $size = $_->[1]-$_->[0];
		my $mid = $_->[0]+floor($size/2);
		$avgpos += $mid;
		$avgcov += $size;
	}
	foreach(@{$hits{$name}->{'d'}}) {
		$avgdiv += $_;
	}
	$avgpos /= $n;
	$avgcov /= $n;
	$avgdiv /= $n;
	print STDERR "\t\t$name\t$n ($frac)\t$avgpos\t$avgcov\t$avgdiv\n";
}

my $Nnone = scalar(@{$repInfo{'none'}});
print STDERR "\tSelecting bg regions...\n";
$c = 0;
foreach(keys %peakInfo) {
	$c++;
	if ($c % 1000 == 0) {
		print STDERR "\t\t$c\n";
	}
	my $pid = $_;
	my $r = $peakInfo{$pid}->{'r'};
	my $rname = $r->{'n'};
	#my $x = {cov=>$cov, n=>$name,coords=>\@coords,div=>$div};
	my $div = $r->{'div'};
	my $cov = $r->{'cov'};
	my $halfSize = floor(($peakInfo{$pid}->{'e'} - $peakInfo{$pid}->{'s'})/2);


	my $startOffset = 0;
	my $endOffset = 2*$halfSize;
	if ($rname ne 'none') {
		$startOffset = $r->{'coords'}->[0];
		$endOffset = $r->{'coords'}->[1];
	}
	#print STDERR "\t$pid\t$rname\t($startOffset,$endOffset)\t$cov\t$div\n";

	#my $r = {s=>$s,e=>$e,ogs=>$s,oge=>$e,c=>$chr,d=>'+',div=>0,
	#	len=>$avgSize,oglen=>$avgSize,frac=>1,n=>$name,mask=>0,sel=>0,score=>0};

	if ($rname eq 'none') {
		for (my $i=0;$i<$factor;$i++) {
			my $ok = 0;
			my $cc = 0;
			my $index = 0;
			while ($ok == 0) {
				$cc++;
				$index = floor($Nnone*rand());
				if ($repInfo{$rname}->[$index]->{'mask'} || $repInfo{$rname}->[$index]->{'sel'}) {
					last if ($cc > 100);
					next;
				}
				$ok = 1;
			}
			next if ($ok ==0);
			my $name = $repInfo{$rname}->[$index]->{'n'};
			my $chr = $repInfo{$rname}->[$index]->{'c'};
			my $s = $repInfo{$rname}->[$index]->{'s'};
			my $e = $repInfo{$rname}->[$index]->{'e'};
			my $d = $repInfo{$rname}->[$index]->{'d'};
			print "$name\t$chr\t$s\t$e\t$d\n";
			$repInfo{$rname}->[$index]->{'sel'} = 1;
		}
		next;
	}


	my $pcov = $endOffset - $startOffset;
	$pcov = 1 if ($pcov < 1);
	my $mid = floor(($endOffset+$startOffset)/2);

	if ($method eq 'sample') {
		my $Nrepeats = scalar(@{$repInfo{$rname}});
		for (my $i=0;$i<$factor;$i++) {
			my $ok = 0;
			my $cc = 0;
			my $index = 0;
			while ($ok == 0) {
				$cc++;
				last if ($cc > 1000*$factor);
				$index = floor($Nrepeats*rand());
				my $r = $repInfo{$rname}->[$index];
				if ($r->{'mask'} || $r->{'sel'}) {
					next;
				}
				my $rstartOffset = $r->{'s'}-$r->{'ogs'};
				my $rendOffset = $r->{'e'}-$r->{'ogs'};
				if ($r->{'d'} eq '-') {
					$rstartOffset = $r->{'oge'}-$r->{'e'};
					$rendOffset = $r->{'oge'}-$r->{'s'};
				}
				if ($rstartOffset < $startOffset) {
					$rstartOffset = $startOffset;
				}
				if ($rendOffset > $endOffset) {
					$rendOffset = $endOffset;
				}
				my $overlap = $rendOffset - $rstartOffset;
				$overlap = 0 if ($overlap < 0);
				my $frac = $overlap / $pcov;
				if ($frac > $overlapThreshold) {
					$ok = 1;
					last;
				}
			}
			next if ($ok ==0);
			my $r = $repInfo{$rname}->[$index];
			$r->{'sel'} = 1;

			my $name = $r->{'n'};
			my $chr = $r->{'c'};
			my $d = $r->{'d'};
			my $m = $r->{'ogs'}+$mid;
			if ($d eq '-') {
				$m = $r->{'oge'}-$mid;
			}
			my $s = $m-$halfSize;
			my $e = $m+$halfSize;
			print "$name\t$chr\t$s\t$e\t$d\n";
		}
		next;
	}

	foreach(@{$repInfo{$rname}}) {
		if ($_->{'mask'} || $_->{'sel'}) {
			$_->{'score'} = -1;
			next;
		}
		if ($_->{'n'} eq 'none') {
			$_->{'score'} = 1;
			next;
		}
		my $rstartOffset = $_->{'s'}-$_->{'ogs'};
		my $rendOffset = $_->{'e'}-$_->{'ogs'};
		if ($_->{'d'} eq '-') {
			$rstartOffset = $_->{'oge'}-$_->{'e'};
			$rendOffset = $_->{'oge'}-$_->{'s'};
		}
		if ($rstartOffset < $startOffset) {
			$rstartOffset = $startOffset;
		}
		if ($rendOffset > $endOffset) {
			$rendOffset = $endOffset;
		}
		my $overlap = $rendOffset - $rstartOffset;
		$overlap = 0 if ($overlap < 0);
		my $frac = $overlap / $pcov;
		$_->{'score'} = $frac;
	}
	my @r = sort {$b->{'score'} <=> $a->{'score'} || 
					$b->{'rand'} <=> $a->{'rand'} } @{$repInfo{$rname}};
	for (my $i=0;$i<$factor;$i++) {
		if ($i >= @r) {
			print STDERR "\t\t\t! ran out of bg repeats for $rname ($pid)\n";
			last;
		}
		my $name = $r[$i]->{'n'};
		my $chr = $r[$i]->{'c'};
		my $d = $r[$i]->{'d'};
		my $m = $r[$i]->{'ogs'}+$mid;
		if ($d eq '-') {
			$m = $r[$i]->{'oge'}-$mid;
		}
		my $s = $m-$halfSize;
		my $e = $m+$halfSize;
		print "$name\t$chr\t$s\t$e\t$d\n";
		$r[$i]->{'sel'} = 1;
	}
}
