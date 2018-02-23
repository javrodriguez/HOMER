#!/usr/bin/env perl
use warnings;



my $foldThresh = 2;
my $fdrThresh = 0.05;
my $peakFoldInput = 2;
my $peakFdrInput = 0.001;
my $style = "";

sub printCMD {
	print STDERR "\n\tUsage: getDifferentialPeaksReplicates.pl [options] -t <IP tagdir1> [IP tagdir2] ...\n";
	print STDERR "\t                                         -b <background tagdir1> [background tagdir2] ...\n";
	print STDERR "\t                                         -i <Input tagdir1> [Input tagdir2] ...\n";
	print STDERR "\t\tNote: if input is provided, peaks will be called.\n";
	print STDERR "\n\tOptions:\n";
	#print STDERR "\t\t-F <#> (fold enrichment over bg, default: $foldThresh)\n";
	#print STDERR "\t\t-fdr <#> (FDR over input, default: $fdrThresh)\n";
	print STDERR "\t\t-f <#> (fold enrichment over bg, default: $foldThresh)\n";
	print STDERR "\t\t-q <#> (FDR over bg, default: $fdrThresh)\n";
	print STDERR "\t\t-fdr <#>, -F <#>, -L <#> (parameters for findPeaks)\n";
	print STDERR "\t\t-genome <genome version> (genome version to use for annotation)\n";
	print STDERR "\t\t-DESeq2 | -DESeq | -edgeR (differential stats algorithm, default: DESeq2)\n";
	print STDERR "\t\t-balanced (normalize signal across peaks, default: normalize to read totals)\n";
	print STDERR "\n\tPeak finding directives:\n";
	print STDERR "\t\t-style <factor|histone|tss> (findPeaks style to use for finding peaks)\n";
	print STDERR "\t\t-use <peaks.txt|regions.txt|tss.txt> (use existing peaks in tag directories)\n";
	print STDERR "\t\t-p <peak file> (use specific peak file instead of tagDir/peaks.txt or finding new one)\n";
	print STDERR "\t\tOther options will be passed to findPeaks\n";
	print STDERR "\n";
	exit;
	
}


my @targets = ();
my @background = ();
my @inputs = ();
my $findPeaksOpts = "";
my $use = "";
my $givenPeakFile = '';
my $norm2total = "-norm2total";
my $diffAlg = "-DESeq2";
my $genome = 'none';
my $annOptions = "";

for (my $i=0;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-t') {
		$i++;
		while ($i < @ARGV) {
			if ($ARGV[$i] =~ /^-/) {
				$i--;
				last;
			}
			push(@targets, $ARGV[$i++]);
		}
	} elsif ($ARGV[$i] eq '-i') {
		$i++;
		while ($i < @ARGV) {
			if ($ARGV[$i] =~ /^-/) {
				$i--;
				last;
			}
			push(@inputs, $ARGV[$i++]);
		}
	} elsif ($ARGV[$i] eq '-b') {
		$i++;
		while ($i < @ARGV) {
			if ($ARGV[$i] =~ /^-/) {
				$i--;
				last;
			}
			push(@background, $ARGV[$i++]);
		}
	} elsif ($ARGV[$i] eq '-f') {
		$foldThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-genome') {
		$genome = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-q') {
		$fdrThresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-use') {
		$use = $ARGV[++$i];
		if ($use eq 'tss.txt') {
			$annOptions .= " -fragLength 1 -strand +";
		}
	} elsif ($ARGV[$i] eq '-p') {
		$givenPeakFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-style') {
		$style = $ARGV[++$i];
		if ($style eq 'tss') {
			$annOptions .= " -fragLength 1 -strand +";
		}
	} elsif ($ARGV[$i] eq '-edgeR') {
		$diffAlg = $ARGV[$i];
	} elsif ($ARGV[$i] eq '-DESeq2') {
		$diffAlg = $ARGV[$i];
	} elsif ($ARGV[$i] eq '-DESeq') {
		$diffAlg = $ARGV[$i];
	} elsif ($ARGV[$i] eq '-balanced') {
		$norm2total = "";
	} elsif ($ARGV[$i] eq '-h' || $ARGV[$i] eq '--help' || $ARGV[$i] eq '--') {
		printCMD();
	} else {
		$findPeaksOpts .= " " . $ARGV[$i];
		#print STDERR "!!! \"$ARGV[$i]\" not recognized\n";
		#printCMD();
	}
}
my $rand = rand();
my %toDelete = ();
if ($diffAlg eq '-edgeR' && $norm2total ne '') {
	print STDERR "!!! Error, -edgeR requires \"-balanced\" to work correctly!!!\n";
	exit;
}
$log2Thresh = log($foldThresh)/log(2);
if (@targets < 1) {
	print STDERR "!!! Error, need at least one target directory!!!\n";
	printCMD();
}
my $targetDirs = '';
my $targetStr = "";
foreach(@targets) {
	$targetDirs .= " \"$_\"";
	$targetStr .= " target";
}
my $inputDirs = '';
my $inputStr = "";
foreach(@inputs) {
	$inputDirs .= " \"$_\"";
	$inputStr .= " input";
}
my $bgDirs = '';
my $bgStr = "";
foreach(@background) {
	$bgDirs .= " \"$_\"";
	$bgStr .= " bg";
}

if ($givenPeakFile eq '') {
	$peakFile = $rand . ".peaks";
	$toDelete{$peakFile}=1;
} else {
	$peakFile = $givenPeakFile;
}

if ($findPeaksOpts ne '') {
	print STDERR "\tUsing the following extra parameters for findPeaks:\n\t\t$findPeaksOpts\n";
}

print STDERR "\tStep1: Defining Putative Peak Set\n";
if ($use eq '' && $givenPeakFile eq '') {
	print STDERR "\t\tFinding peaks in merged meta-experiment from target tag directories\n";
	my $targetDir = $rand . ".targetTagDir";
	`makeTagDirectory \"$targetDir\" -d $targetDirs`;
	$toDelete{$targetDir}=1;
	my $inputDir = $rand . ".inputTagDir";
	my $cmd = "findPeaks \"$targetDir\" $findPeaksOpts";
	if ($style eq '') {
		$style = 'factor';
		print STDERR "\tUsing -style $style...\n";
	}
	$cmd .= " -style $style";
	if (@inputs > 0) {
		`makeTagDirectory \"$inputDir\" -d $inputDirs`;
		$toDelete{$inputDir}=1;
		$cmd .= " -i \"$inputDir\"";
	}
	#print STDERR "`$cmd > \"$peakFile\"`\n";
	`$cmd > \"$peakFile\"`;
} elsif ($givenPeakFile eq '' && $use ne '') {
	my $files = "";
	my @allDirs = ();
	push(@allDirs, @targets, @inputs, @background);
	print STDERR "\t\tUsing existing peak files for features:\n";
	foreach(@allDirs) {
		my $p = $_ . "/" . $use;
		if (-e $p) {
			print STDERR "\t\t\t$p\n";
			$files .= " \"$p\"";
		}
	}
	`mergePeaks $files > \"$peakFile\"`;
}



$rawFile= $rand . ".raw.txt";
$diffFile= $rand . ".diff.txt";
$upFile = $rand . ".Up_target_vs_bg.txt";
$downFile = $rand . ".Down_target_vs_bg.txt";
if ($bgDirs eq '') {
	$bgDirs = $inputDirs;
	$bgStr = $inputStr;
	@background = @inputs;
	$upFile = $rand . ".Up_target_vs_input.txt";
	$downFile = $rand . ".Down_target_vs_input.txt";
}

print STDERR "\n\tStep2: Quantifying reads across target/background/input tag directories\n\n";
#print STDERR "`annotatePeaks.pl \"$peakFile\" none -d $bgDirs $targetDirs -raw > \"$rawFile\"`;\n";
`annotatePeaks.pl \"$peakFile\" $genome -d $bgDirs $targetDirs -raw $annOptions > \"$rawFile\"`;
#print STDERR "`getDiffExpression.pl \"$rawFile\" $bgStr $targetStr $norm2total $diffAlg -fdr $fdrThresh -log2fold $log2Thresh -export $rand > $diffFile`;\n";
#

print STDERR "\n\tStep3: Calling R for differential enrichment statistics ($diffAlg)\n\n";
`getDiffExpression.pl \"$rawFile\" $bgStr $targetStr $norm2total $diffAlg -fdr $fdrThresh -log2fold $log2Thresh -export $rand > $diffFile`;
$toDelete{$rawFile}=1;
$toDelete{$diffFile}=1;
$toDelete{$upFile}=1;
$toDelete{$downFile}=1;
open IN, $upFile;
while (<IN>) {
	print $_;
}
close IN;

foreach(keys %toDelete) {
	next if ($_ eq '/');
	#print STDERR "\trm -r \"$_\"\n";
	`rm -r "$_"`;
}


