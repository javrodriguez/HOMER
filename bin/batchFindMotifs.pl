#!/usr/bin/env perl
use warnings;



if (@ARGV < 3) {
	print STDERR "\n\tUsage: batchFindMotifs.pl [promoter set] [options...] -f list1.txt list2.txt ...\n";
	print STDERR "\n\t\t-cpu (# of concurrent jobs, -p controls CPUs used by each findMotifsGenome.pl instance)\n";
    print STDERR "\n";
	print STDERR "\n";
	exit;
}

my $promoterSet = $ARGV[0];
my $options = "";
my @files = ();
my $maxCPUs = 1;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-f') {
		for (my $j=$i+1;$j<@ARGV;$j++) {
			push(@files, $ARGV[$j]);
		}
		last;
	} elsif ($ARGV[$i] eq '-cpu') {
		$maxCPUs = $ARGV[++$i];
	} else {
		$options .= " " . $ARGV[$i];
	}
}

my $cpus = 0;
my @pids = ();
foreach(@files) {
	my $dir = $_;
	$dir =~ s/\.txt$//;
	$dir =~ s/\.bed$//;
	$dir = "Motifs-$dir";

	my $pid = fork();
	$cpus++;
	if ($pid == 0) {
		print STDERR "findMotifs.pl \"$_\" $promoterSet \"$dir\" $options\n";
		`findMotifs.pl "$_" $promoterSet "$dir" $options`;
		exit(0);
	}
	push(@pids, $pid);
    if ($cpus >= $maxCPUs) {
        my $id = wait();
        print "\t$id finished\n";
        $cpus--;
    }
}

my $id = 0;
while ($id >= 0) {
    $id = wait();
    if ($id == -1) {
        print STDERR "\tALL FINISHED!!!\n";
    } else {
        print STDERR "\t$id finished\n";
    }
}
print STDERR "\n";
