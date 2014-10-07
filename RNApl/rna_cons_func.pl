# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;

use List::Util qw(sum);

# 
sub getTargetSiteConservationFeatures(\@\@$)
{
	my @targetCandidate = @{shift()};
	my @phastConsScores = @{shift()};
	my $flankLength = shift() || 70;
	my $utrLength = scalar(@phastConsScores);
	my $startPos = $targetCandidate[1] - 1;
	my $endPos = $targetCandidate[2] - 1;
	my $upStart = $startPos - $flankLength;
	my $downEnd = $endPos + $flankLength;
	my $seedScore = 0;
	my $siteScore = 0;
	my $flankScore = 0;
	
	if ($upStart < 0)
	{
		$upStart = 0;
	}
	
	if ($downEnd >= $utrLength)
	{
		$downEnd = $utrLength - 1;
	}
	
	#print STDERR "$utrLength $startPos $endPos\n";
	
	if ($endPos < $utrLength)
	{
		# Seed score
		my $startSeed = $endPos - 7;
		my @seedScores = @phastConsScores[$startSeed .. $endPos];
		$seedScore = sum(@seedScores) / @seedScores;
		
		#print STDERR "$seedScore\n";
		
		# Target site score
		my @siteScores = @phastConsScores[$startPos .. $endPos];
		$siteScore = sum(@siteScores) / @siteScores;
		
		# Flank score
		my @flankScores = @phastConsScores[$upStart .. ($startPos - 1)];
		my @downScores = @phastConsScores[($endPos + 1) .. $downEnd];
		push(@flankScores, @downScores);
		$flankScore = sum(@flankScores) / @flankScores;
	}
	
	return($seedScore, $siteScore, $flankScore);
}

1;
