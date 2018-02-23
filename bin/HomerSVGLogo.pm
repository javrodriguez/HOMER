#!/usr/bin/env perl
use warnings;


package HomerSVGLogo;

my $defWidth=25;
my $defHeight=50;
my $defFont = "Arial";

########################################
# motif2svg($motif,$bigFlag,$rnaFlag)
#
# $motif is a hash containg the following:
#   $motif->{'matrix'} : 2D array (length x nucleotide) with frequencies
# $bitFlag = 1 (for information content scaling) or 0 (frequency scaling)
# $rnaFlag = 1 (for U instead of T) or 0 (DNA letters)
# 
# Returns a string that contains the SVG code (either put directly in HTML or SVG file
########################################

sub motif2svg {
	my ($motif,$bitFlag,$rnaFlag) = @_;
	my $str = "";

	# here we set all of the specific settings to make the logos look OK with Arial font...
	my @alpha =();
	$alpha[0]='A';
	$alpha[1]='C';
	$alpha[2]='G';
	$alpha[3]='T';
	my @color =();
	$color[0]='#00BB00';
	$color[1]='#0000EE';
	$color[2]='#F9A500';
	$color[3]='#DD0000';
	my @xwidth =();
	$xwidth[0]=0.53;
	$xwidth[1]=0.59;
	$xwidth[2]=0.55;
	$xwidth[3]=0.65;
	my @xoffset=();
	$xoffset[0]=2;
	$xoffset[1]=0;
	$xoffset[2]=0;
	$xoffset[3]=1;

	if ($rnaFlag == 1) {
		$alpha[3]='U';
		$xoffset[3]=-1.5;
	}

	my $motifLength = scalar(@{$motif->{'matrix'}});

	my $finalLength = $defWidth*$motifLength+5;

	$str .= "<svg width=\"$finalLength\" height=\"$defHeight\">\n";
	$str .= " <g font-family=\"$defFont\" font-weight=\"bold\" font-size=\"66.5\">\n";

	my @matrix = @{$motif->{'matrix'}};

	if ($bitFlag) {
		my @bitMatrix = ();
		my $ecorr = 0;

		for (my $i=0;$i<@matrix;$i++) {
			my @bits = ();
			my $H = 0;
			for (my $j=0;$j<@{$matrix[$i]};$j++) {
				my $v = $matrix[$i][$j];
				$v = 0.001 if ($v < 0.001);
				$H += $v * log($v)/log(2);
			}
			$H *= -1;
			for (my $j=0;$j<@{$matrix[$i]};$j++) {
				$bits[$j] = $matrix[$i][$j] * (2 - ($H+$ecorr));

				#scale bit content by 2 so that it can be plotted the same way as the probabilities
				$bits[$j] /= 2;
			}
			push(@bitMatrix, \@bits);
		}
		@matrix = @bitMatrix;
	}
		


	for (my $i=0;$i<@matrix;$i++) {

		my @order = (0,1,2,3);
		@order = sort {$matrix[$i][$a] <=> $matrix[$i][$b]} @order;

		my $currX = $i*$defWidth;
		my $currY = $defHeight-1;

		foreach(@order) {
			my $prob = $matrix[$i][$_];
			my $nuc = $alpha[$_];
			my $color = $color[$_];
			my $xscale = $xwidth[$_];

			my $yscale = $prob;
			$xscale = $xscale;

			my $x = $currX + $xoffset[$_];
			my $y = $currY;

			$str .= "  <text fill=\"$color\" x=\"0\" y=\"0\" ";
			$str .= " transform=\"matrix($xscale,0,0,$yscale,$x,$y)\">";
			$str .= "$nuc</text>\n";

			$currY -= $prob*$defHeight;

		}
		$str .= "  \n";
	}

	$str .= " </g>\n";
	$str .= "</svg>\n";

	return $str;
}
1;

