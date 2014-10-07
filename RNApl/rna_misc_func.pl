# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;
use List::Util qw(min);

# Convert RNA sequence to A, C, G, and U
sub formatRnaSeq($)
{
	my $rnaSeq = shift();
	
	$rnaSeq =~ s/a/A/g;
	$rnaSeq =~ s/c/C/g;
	$rnaSeq =~ s/g/G/g;
	$rnaSeq =~ s/t/U/g;
	$rnaSeq =~ s/T/U/g;
	$rnaSeq =~ s/u/U/g;
	
	return $rnaSeq;
}

# Return distance of target site to closest 3' UTR end
sub getTargetSiteEndDistance($\@)
{
	my $utrSeq = shift();
	my @targetCandidate = @{shift()};
	
	my @dist = ($targetCandidate[1], 
	            length($utrSeq) - $targetCandidate[2]);
	
	return(min(@dist) / length($utrSeq));
}

1;
