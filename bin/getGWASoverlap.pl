#!/usr/bin/env perl
use warnings;

use POSIX;
use Statistics;
use File::Spec;
use File::Basename;

my $dist = 1000;
my $minSNPs =10;
my $maxCPUs = 1;
my $gregorPopulation = "EUR";
my $LDthresh = 0.8;

sub printCMD {
	print STDERR "\n\tgetGWASoverlap.pl <gwas catolog file> -p <peak file1> [peak file2] ... [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-d <#> (Overlap distance, default: $dist)\n";
	print STDERR "\t\t-min <#> (minimum number of significant SNPs to consider study, default: $minSNPs)\n";
	print STDERR "\t\t-cpu <#> (number of threads to use, default: $maxCPUs)\n";
	#print STDERR "\t\t-empirical <#> (calculate FDR q-values empirically, <#> number of randomizations)\n";
	print STDERR "\t\t-GREGOR <path-to-GREGOR> <path-to-Reference LD info> (perform enrichment test with GREGOR)\n";
	print STDERR "\t\t\t-LD <#> (LD threshold for 'buddy SNPs', default 0.8, must be greater than 0.7)\n";
	print STDERR "\t\t-pubmedID <#> (only analyze this study)\n";
	print STDERR "\t\t-snpOut <output file> (output overlapping snps & buddies as a BED file)\n";
	print STDERR "\n\tThe gwas catalog file can be downloaded from UCSC annotation database:\n";
	print STDERR "\t\ti.e.: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gwasCatalog.txt.gz\n";
	print STDERR "\n";
	exit;
}



if (@ARGV < 1) {
	printCMD();
}


my $rand = rand();
my $gregorFlag = 0;

my $gwasFile = $ARGV[0];
my @peakFiles = ();
my $numEmpirical = 0;
my $pubmedID = '';
my $snpOutFile = '';

for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-d') {
		$dist = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-min') {
		$minSNPs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-LD') {
		$LDthresh = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-snpOut') {
		$snpOutFile = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-pubmedID') {
		$pubmedID = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-empirical') {
		$numEmpirical = $ARGV[++$i];
		print STDERR "!!! Error - '-empirical' is disabled in this version\n";
		exit;
	} elsif ($ARGV[$i] eq '-GREGOR' || $ARGV[$i] eq '-gregor') {
		$gregorFlag = 1;
		$gregorExe = $ARGV[++$i];
		$gregorRef = $ARGV[++$i];
		print STDERR "\tWill use GREGOR:\n";
		print STDERR "\t\tprogram location: $gregorExe\n";
		print STDERR "\t\tLD info location: $gregorRef\n";
	} elsif ($ARGV[$i] eq '-p') {
		$i++;
		while ($i<@ARGV && $ARGV[$i] !~ /^-/) {
			push(@peakFiles, $ARGV[$i++]);
		}
		$i--;
	} else {
		print STDERR "Couldn't recognize $ARGV[$i]\n";
		printCMD();
	}
}

if ($gregorFlag) {
	$maxCPUs = ceil($maxCPUs/10.0);
	print STDERR "\tAdjusting simultaneous processes to $maxCPUs (GREGOR uses ~10 cpus per process)\n";
}
my %limits = ();


print STDERR "\n\tMinimum Risk SNPs per study: $minSNPs\n";
print STDERR "\tMaximum Distance between Peaks and risk SNPs: $dist\n\n";

print STDERR "\n\tAnalyzing peak files:\n";
my @fdrPvalues = ();

for (my $i=0;$i<@peakFiles;$i++) {
	print STDERR "\t\t$peakFiles[$i]\n";
	if ($numEmpirical > 0) {
		my $tmpFile = $rand . ".tmp";
		`bed2pos.pl "$peakFiles[$i]" -check > "$tmpFile"`;
		open IN, $tmpFile;
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			next if ($line[0] =~ /^#/);
			next if (@line < 5);
			if (!exists($limits{$line[1]})) {
				$limits{$line[1]} =0;
			}
			$limits{$line[1]} = $line[3] if ($limits{$line[1]} < $line[3]);
		}
		close IN;
		`rm $tmpFile`;
		my @a = ();
		push(@fdrPvalues, \@a);
	}
}
print STDERR "\n";

my @bedFiles = ();
my $gregorBedInputIndexFile = $rand . ".gregor.bed.index";
my $absBedIndexFile = "";
my $absGregorRef = "";
if ($gregorFlag) {
	open GREGINDEX, ">$gregorBedInputIndexFile";
	for (my $i=0;$i<@peakFiles;$i++) {
		my $tmpFile = $rand . ".tmp";
		`bed2pos.pl "$peakFiles[$i]" -check > "$tmpFile"`;
		my $bedFileName= $rand . "." . $peakFiles[$i] . ".bed";
		`pos2bed.pl "$tmpFile" > "$bedFileName"`;
		push(@bedFiles, $bedFileName);
		my $absFile = File::Spec->rel2abs($bedFileName);
		print GREGINDEX "$absFile\n";
	}
	close GREGINDEX;
	$absBedIndexFile = File::Spec->rel2abs($gregorBedInputIndexFile);
	$absGregorRef = File::Spec->rel2abs($gregorRef);
}
my $nbed = scalar(@bedFiles);


#print out the header
print "Study";
print "\tTotal risk SNPs";
for (my $i=0;$i<@peakFiles;$i++) {
	print "\t$peakFiles[$i] Overlap";
	print "\t$peakFiles[$i] p-value";
	if ($numEmpirical > 0) {
		print "\t$peakFiles[$i] q-value FDR/Empirical(n=$numEmpirical)";
	} else {
		print "\t$peakFiles[$i] q-value FDR/Benjamini";
	}
	print "\t$peakFiles[$i] SNPs";
}
print "\n";

my $limitTotal = 0;
if ($numEmpirical > 0) {
	#print STDERR "\tChromosome Limits:\n";
	foreach(keys %limits) {
		#print STDERR "\t\t$_\t$limits{$_}\n";
		$limitTotal += $limits{$_};
	}
	#print STDERR "\tApprox Size: $limitTotal\n\n";
}

my $gwas = parseGWAScatalog($gwasFile, $pubmedID);


@results = ();
@studies = ();

my $numStudy = scalar(keys %$gwas);
my $count2 = 0;
my $cpus=0;
my %tmpIDs = ();
foreach(keys %$gwas) {
	$count2++;
	#last if ($count2 > 40);
	my $study = $_;
	my @snps = keys %{$gwas->{$study}};
	my $N = scalar(@snps);
	next if (@snps < $minSNPs);

	my $tmpID = $rand . "-$count2";
	my $pid = fork();
	if ($pid != 0) {
		#parent process
		#print STDERR "$pid\t$tmpID\n";
		$tmpIDs{$pid} = $tmpID;
		$cpus++;
		if ($cpus >= $maxCPUs) {
			my $id = wait();
			my $tid = $tmpIDs{$id};
			#print STDERR "\tProcessing tmpID: $tid\t process ID: $id\n";
			processResults($tid);
			$cpus--;
		}
		next;
	}


	print STDERR "\tAnalyzing Study $count2 of $numStudy ($study, $N total risk SNPs)\n";

	my @peakResults = ();

	my $resultsFile = $tmpID . ".results";
	open RESULTS, ">$resultsFile";

	my $numDsnps = $N;
	print RESULTS "$study\t$N\t$tmpID\n";
	#push(@studies, "$study\t$N");


	if ($gregorFlag) {
		my $snpIndexFile = "$tmpID" . ".snps.txt";
		open OUT, ">$snpIndexFile";
		foreach(@snps) {
			my $allele = $_;
			$allele =~ s/\-.+$//;
			#my $p = $gwas->{$study}->{$allele};
			print OUT "$allele\n";
		}
		close OUT;
		my $absSnpIndexFile = File::Spec->rel2abs($snpIndexFile);
		my $outputDir = $tmpID . "-outDir";
	   	$absOutputDir = File::Spec->rel2abs($outputDir);
		
		my $confFile = $tmpID . ".conf";
		open OUT, ">$confFile";
		print OUT "#############################################################################\n";
		print OUT "# CHIPSEQ ENRICHMENT CONFIGURATION FILE\n";
		print OUT "# This configuration file contains run-time configuration of\n";
		print OUT "# CHIP_SEQ ENRICHMENT\n";
		print OUT "#############################################################################\n";
		print OUT "# KEY ELEMENTS TO CONFIGURE : NEED TO MODIFY\n";
		print OUT "#############################################################################\n";
		print OUT "INDEX_SNP_FILE = $absSnpIndexFile\n";
		print OUT "BED_FILE_INDEX = $absBedIndexFile\n";
		print OUT "REF_DIR = $absGregorRef\n";
		print OUT "R2THRESHOLD = $LDthresh ## must be greater than 0.7\n";
		print OUT "LDWINDOWSIZE = 1000000 ## must be less than 1MB; these two values define LD buddies\n";
		print OUT "OUT_DIR = $absOutputDir\n";
		print OUT "MIN_NEIGHBOR_NUM = 500 ## define the size of neighborhood\n";
		print OUT "BEDFILE_IS_SORTED = false  ## false, if the bed files are not sorted\n";
		print OUT "POPULATION = $gregorPopulation  ## define the population, you can specify EUR, AFR, AMR or ASN\n";
		print OUT "TOPNBEDFILES = 2000\n";
		print OUT "JOBNUMBER = 10\n";
		print OUT "#############################################################################\n";
		print OUT "#############################################################################\n";
		print OUT "BATCHTYPE = local ##  run jobs on local machine\n";
		close OUT;

		#print STDERR "\tRunning GREGOR\n";
		#print STDERR "perl $gregorExe --conf $confFile\n";
		`perl $gregorExe --conf $confFile`;
		#print STDERR "\tGREGOR finished\n";

		open IN, "$outputDir/StatisticSummaryFile.txt";
		my $c = 0;
		my %peaks = ();
		while (<IN>) {
			$c++;
			next if ($c==1);
			chomp;
			s/\r//g;
			my @line = split /\t/;
			my $bedname = "$line[0]";
			my $oexp = "$line[1] ($line[2])";
			my $pvalue = $line[3];
			$peaks{$bedname} = {oexp=>$oexp,p=>$pvalue};
		}
		close IN;

		`ls -d $outputDir/index.SNP.and.LD.for.top.*/ > $outputDir/.dir`;
		open IN, "$outputDir/.dir";
		my $dd = $outputDir . "/index.SNP.and.LD.for.top." . $nbed . ".bed/";
		while (<IN>) {
			chomp;
			$dd = $_;
			last;
		}
		close IN;

		foreach(@bedFiles) {
			my $bedname = $_;
			my $N = 0;
			my $pvalue = 1;
			my $snpStr = '';
			my $fdr = 1;
			my $pvalueStr = '';
			if (exists($peaks{$bedname})) {
				$N = $peaks{$bedname}->{'oexp'};
				$pvalue = $peaks{$bedname}->{'p'};
			}
			my $overlapSNPsFile = $dd . "index.SNP.and.LD.in.bed." . $bedname . ".txt";
			if (-e $overlapSNPsFile) {
				print STDERR "Cool - found overlap file.\n";
				open IN, $overlapSNPsFile;
				while (<IN>) {
					chomp;
					s/\r//g;
					my @line = split /\t/;
					next if (@line < 4);
					next if ($line[2] eq '-1');
					next if ($line[1] eq 'LD');
					$snpStr .= "chr" . $line[1] . ",";
				}
				close IN;
			} else {
				print STDERR "Uh... $overlapSNPsFile is not right. Possible:\n";
				print STDERR "Listing of output Directory:\n";
				`ls $outputDir/`;
				`ls $outputDir/index.SNP.and.LD.for.top*`;
			}

			print RESULTS "$bedname\t$N\t$pvalue\t$snpStr\t$fdr\t$pvalueStr\n";
		}
		close RESULTS;

		#print STDERR "sleeping...\n";
		#`sleep 10000`;

		#`rm "$snpIndexFile"`;
		#`rm -r "$outputDir"` if ($outputDir ne '' || $outputDir ne '/');

		exit(0);
	}


	# HOMER simplified enrichment scheme
	my $tmpFile = $tmpID . ".tmp";
	open OUT, ">$tmpFile";
	foreach(@snps) {
		my $allele = $_;
		my $p = $gwas->{$study}->{$allele};
		print OUT "$p\n";
	}
	close OUT;
	my $tmpFile2 = $tmpID . ".2.tmp";

	for (my $i=0;$i<@peakFiles;$i++) {
		`cp "$peakFiles[$i]" "$tmpFile2"`;

		`mergePeaks "$tmpFile" "$tmpFile2" -matrix $tmpID  -d $dist > /dev/null 2> /dev/null`;

		my $logp = 0;

		my $c = 0;
		open IN, "$tmpID.logPvalue.matrix.txt";
		while (<IN>) {
			$c++;
			next if ($c < 2);
			s/\r//g;
			my @line = split /\t/;
			$logp = $line[2];
			last;
		}
		close IN;
		my $pvalue = 1;
		if ($logp < 0.1) {
			$pvalue = exp($logp);
		} else {
			$pvalue = 1-exp(-1*$logp);
		}
		`rm $tmpID.*matrix.txt`;

		`mergePeaks "$tmpFile" "$tmpFile2" -cobound 1 -prefix $tmpID -d $dist > /dev/null 2> /dev/null`;

		my @bound = ();
		my $snpStr = "";
		open IN, "$tmpID.coBoundBy1.txt";	
		while (<IN>) {
			chomp;
			next if (/^#/);
			my @line = split /\t/;
			push(@bound, $line[0]);
			$snpStr .= "$line[0]" . ",";
		}
		close IN;
		my $NN = scalar(@bound);	
		`rm $tmpID.coBoundBy*`;

		my $fdr = 1;
		my $pvalueStr = '';

		if ($numEmpirical > 0) {
			my $fileStr = "";
			my @chrs = keys %limits;
			for (my $j=0;$j<$numEmpirical;$j++) {
				my @rands = ();
				for (my $k=0;$k<$numDsnps;$k++) {
					push(@rands, floor($limitTotal*rand()));
					#print STDERR "$j $k $rands[$k]\n";
				}
				@rands = sort {$a <=> $b} @rands;
				my $f = "$tmpID.empirical.$j.tmp";
				$fileStr .= " $f";
				open OUT, ">$f";
				my $chrIndex = 0;
				my $totalIndex = 0;
				foreach(@rands) {
					my $r = $_;	
					while ($totalIndex+$limits{$chrs[$chrIndex]}< $r) {
						$totalIndex += $limits{$chrs[$chrIndex++]};
					}
					my $p = $r-$totalIndex;
					my $p2 = $p+1;
					print OUT "$r\t$chrs[$chrIndex]\t$p\t$p2\t+\n";
				}
				close OUT;
			}
			#print STDERR "`mergePeaks $tmpFile2 -cobound 1 -matrix $tmpID -d $dist $fileStr > /dev/null`;\n";
			`mergePeaks "$tmpFile2" -cobound 1 -matrix $tmpID -d $dist $fileStr -gsize $limitTotal > /dev/null 2> /dev/null`;

			my @pvalues = ();
			open IN, "$tmpID.logPvalue.txt";
			my $c = 0;
			while (<IN>) {
				$c++;
				next if ($c < 2);
				chomp;
				my @line = split /\t/;
				my $p = 1;
				if ($line[1] <= 0.1) {
					$p = exp($line[1]);
				} else {
					$p = 1-exp(-1*$line[1]);
				}
				push(@pvalues, $p);
				$pvalueStr .= "$p|";
				#print STDERR "\t$pvalue\t$p\n";
			}
			close IN;
			$pvalueStr .= ",";
			`rm -f $fileStr`;
			push(@{$fdrPvalues[$i]},\@pvalues);

		}

		#my $d = {n=>$N,p=>$pvalue,snps=>$snpStr,f=>$fdr};	
		#push(@peakResults, $d);
		print RESULTS "$peakFiles[$i]\t$NN\t$pvalue\t$snpStr\t$fdr\t$pvalueStr\n";
	}
	#push(@results, \@peakResults);
	`rm "$tmpFile" "$tmpFile2"`;

	close RESULTS;
	exit(0);
}

$id = 0;
while ($id >= 0) {
	$id = wait();
	last if ($id < 0);
	my $tid = $tmpIDs{$id};
	#print STDERR "\tProcessing tmpID: $tid\t process ID: $id\n";
	processResults($tid);
}


for (my $i=0;$i<@peakFiles;$i++) {
	my @pvalues = ();
	foreach(@results) {
		push(@pvalues, $_->[$i]->{'p'});
	}
	my $fdrs = "";
	if ($numEmpirical < 1) {
		$fdrs = Statistics::benjaminiFDR(\@pvalues);
	} else {
		my $aa = scalar(@pvalues);
		my $bb = scalar(@{$fdrPvalues[$i]});
		$fdrs = Statistics::empiricalFDR(\@pvalues, $fdrPvalues[$i]);
	}
	for (my $j=0;$j<@results;$j++) {
		$results[$j]->[$i]->{'f'} = $fdrs->[$j];
	}
}

if ($snpOutFile ne '') {
	open SNPBED, ">$snpOutFile";
}
for (my $i=0;$i<@results;$i++) {
	print $studies[$i];
	my $z = 0;
	foreach(@{$results[$i]}) {
		print "\t$_->{'n'}\t$_->{'p'}\t$_->{'f'}\t$_->{'snps'}";
		if ($snpOutFile ne '') {
			my @snps = split /\,/,$_->{'snps'};
			#print STDERR "$_->{'snps'} = @snps\n";
			my $sname = $studies[$i];
			my $expName = $peakFiles[$z];
			$sname =~ s/\s/_/g;
			foreach(@snps) {
				my $bname = $sname . "|" . $expName . "|$_";
				my @a = split /\:/;
				next if (@a < 2);
				my $s = $a[1]-1;
				print SNPBED "$a[0]\t$s\t$a[1]\t$bname\t0\t+\n";
			}
		}
		$z++;
	}
	print "\n";
}	
if ($snpOutFile ne '') {
	close SNPBED;
}

if ($gregorFlag) {
	foreach(@bedFiles) {
		`rm "$_"`;
	}
}

`rm -f $rand*`;


sub processResults {
	my ($tid) = @_;
	open IN, "$tid.results" or print STDERR "!!! Could not open results for $tid\n";
	my $c = 0;
	my $study = '';
	my $N = 0;
	my @r = ();
	my $peakIndex = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c == 1) {
			$study = $line[0];
			$N = $line[1];
			next;
		}
		#print RESULTS "$peakFiles[$i]\t$N\t$pvalue\t$snpStr\t$fdr\t$pvalueStr\n";
		my $n = $line[1];
		my $pvalue = $line[2];
		my $snpStr = $line[3];
		my $fdr = $line[4];
		my $d = {n=>$n,p=>$pvalue,snps=>$snpStr,f=>$fdr};	
		push(@r, $d);

		my $epi = '';
		if (@line > 5) {
			$epi = $line[5];
			my @randomizations = split /\,/, $epi;
			if (@randomizations > 1) {
				pop @randomizations;
			}
			foreach(@randomizations) {
				my @ps = split /\|/, $_;
				pop @ps;
				push(@{$fdrPvalues[$peakIndex]},\@ps);
			}
		}
		$peakIndex++;

	}
	close IN;
	push(@studies, "$study\t$N");		
	push(@results, \@r);
	`rm "$tid.results"`;
}



sub parseGWAScatalog {
	my ($file, $pubmedID) = @_;
	my %gwas = ();
	my $format = 'gwas';

	open IN, $file;
	my $c = 0;
	while (<IN>) {
		$c++;
		chomp;
		s/\r//g;
		my @line = split /\t/;
		if ($c == 1) {
			if (@line < 20) {
				print STDERR "\tAssuming SNP file foramt\n";
				$format = 'snps';
			} else {
				print STDERR "\tAssuming file is UCSC GWAS Catalog format...\n";
				$format = 'gwas';
			}
		}

		if ($format eq 'gwas') {
			my $chr = $line[1];
			my $start = $line[2];
			my $end = $line[3];
			my $snp = $line[4];
			my $pubmed = $line[5];
			if ($pubmedID ne '') {
				next if ($pubmed ne $pubmedID);
			}
			my $study = $line[10] . "($pubmed)";
			my $allele = $line[15];
			my $pvalue = $line[17];
			my $position = "$allele\t$chr\t$start\t$end\t0";
		
			if (!exists($gwas{$study})) {
				my %snps = ();
				$gwas{$study} = \%snps;
			}
			$gwas{$study}->{$allele} = $position;
		} elsif ($format eq 'snps') {
			my $study = $file;
			my $allele = $line[0];
			my $position = "$allele" . "\t" . $_;
			if (!exists($gwas{$study})) {
				my %snps = ();
				$gwas{$study} = \%snps;
			}
			$gwas{$study}->{$allele} = $position;
		}
	}
	close IN;
	return \%gwas;
}
