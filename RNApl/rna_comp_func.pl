# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;
use Bio::Perl;
use Bio::Tools::SeqStats;
use Bio::Tools::SeqWords;

# Compute composition features
sub getCompositionFeatures($\@$)
{
	my $rnaSeq = shift();
	my @targetCandidate = @{shift()};
	my $flankLength = shift() || 70;
	
	my $targetSeq = substr($rnaSeq, $targetCandidate[1], 
	                       $targetCandidate[2]-$targetCandidate[1]+1);
	my $upStart = $targetCandidate[1] - $flankLength - 1;
	my $downEnd = $targetCandidate[2] + $flankLength - 1;
	
	#print STDERR "$targetCandidate[1] $targetCandidate[2] $upStart $downEnd\n";
	
	if ($upStart < 0)
	{
		$upStart = 0;
	}
	
	if ($downEnd >= length($rnaSeq))
	{
		$downEnd = length($rnaSeq) - 1;
	}
	
	my $upFlankSeq = substr($rnaSeq, $upStart,
	                        $targetCandidate[1] - $upStart);
	my $downFlankSeq = substr($rnaSeq, $targetCandidate[2] - 1, 
	                          $downEnd - $targetCandidate[2]);
	
	my $len = length($rnaSeq);
	#print STDERR "mRNA length: $len\n";
	#print STDERR "Target site: $targetCandidate[1] $targetCandidate[2]\n";
	#print STDERR "Target: $targetSeq\n";
	#print STDERR "Up flank: $upFlankSeq\n";
	#print STDERR "Down flank: $downFlankSeq\n";
	                          
	# Monomer count
	my $targetSeqObj = new_sequence($targetSeq);
	my $targetMonomers = 
		Bio::Tools::SeqStats->count_monomers($targetSeqObj);
	#foreach my $base (sort keys %$targetMonomers) 
	#{
	#	print "Number of bases of type ", $base, "= ", 
	#	$targetMonomers->{$base},"\n";
	#}
	
	my $upFlankMonomers;
	if (length($upFlankSeq) > 0)
	{
		my $upFlankSeqObj = new_sequence($upFlankSeq);
		$upFlankMonomers = 
			Bio::Tools::SeqStats->count_monomers($upFlankSeqObj);
	}
	
	my $downFlankMonomers;
	if (length($downFlankSeq) > 0)
	{
		my $downFlankSeqObj = new_sequence($downFlankSeq);
		$downFlankMonomers = 
			Bio::Tools::SeqStats->count_monomers($downFlankSeqObj);
	}
	
	# Dimer count
	my $targetDimers = 
		Bio::Tools::SeqWords->count_words($targetSeqObj, 2);
	#foreach my $base (sort keys %$targetDimers) 
	#{
	#	print "Number of bases of type ", $base, "= ", 
	#	$targetDimers->{$base},"\n";
	#}
	my $upFlankDimers;
	my $downFlankDimers;
	
	if (length($upFlankSeq) > 1)
	{
		my $upFlankSeqObj = new_sequence($upFlankSeq);
		$upFlankDimers = 
			Bio::Tools::SeqWords->count_words($upFlankSeqObj, 2);
	}
	
	if (length($downFlankSeq) > 1)
	{
		my $downFlankSeqObj = new_sequence($downFlankSeq);
		$downFlankDimers = 
			Bio::Tools::SeqWords->count_words($downFlankSeqObj, 2);
	}
	
	return ($targetMonomers, $upFlankMonomers, $downFlankMonomers,
	        $targetDimers, $upFlankDimers, $downFlankDimers);
}

sub getMonomerList()
{
	return(["A", "C", "G", "U"]);
}

sub getDimerList()
{
	my @alpha = @{getMonomerList()};
	my $n = @alpha;
	my @list = ();
	my $count = 0;
	
	for (my $i = 0; $i < $n; $i++)
	{
		for (my $j = 0; $j < $n; $j++)
		{
			$list[$count++] = $alpha[$i] . $alpha[$j];
		}
	}
	
	return(\@list);
}

sub monomerHash2Array(\%)
{
	my $ref = shift();
	my %hash;
	
	if (defined($ref))
	{
		%hash = %{$ref};
	}
	
	my @alpha = @{getMonomerList()};
	my @array = ();
	my $i = 0;
	my $total = 0;
	
	foreach my $key (@alpha)
	{
		$array[$i] = 0;
		
		if (defined($hash{$key}))
		{
			#$array[$i++] = $hash{$key};
			$array[$i] = $hash{$key};
			$total += $hash{$key};
		}
		#else
		#{
		#	$array[$i++] = 0;
		#}
		
		$i++;
	}
	
	if ($total > 0)
	{
		for(my $j = 0; $j < scalar(@array); $j++)
		{
			$array[$j] /= $total;
		}
	}
	
	return(\@array);
}

sub dimerHash2Array(\%)
{
	my $ref = shift();
	my %hash;
	
	if (defined($ref))
	{
		%hash = %{$ref};
	}
	
	my @alpha = @{getDimerList()};
	my @array = ();
	my $i = 0;
	my $total = 0;
	
	foreach my $key (@alpha)
	{
		$array[$i] = 0;
		
		if (defined($hash{$key}))
		{
			#$array[$i++] = $hash{$key};
			$array[$i] = $hash{$key};
			$total += $hash{$key};
		}
		#else
		#{
		#	$array[$i++] = 0;
		#}
		
		$i++;
	}
	
	if ($total > 0)
	{
		for(my $j = 0; $j < scalar(@array); $j++)
		{
			$array[$j] /= $total;
		}
	}
	
	return(\@array);
}

# AU in seed flanks based on Grimson et al. 2007
sub getSeedFlankAUContent($\@$$)
{
	my $rnaSeq = shift();
	my @targetCandidate = @{shift()};
	my $seedLength = shift() || 8;
	my $flankLength = shift() || 30;
	
	my $upStart = $targetCandidate[2] - $seedLength - $flankLength - 1;
	my $upEnd = $targetCandidate[2] - $seedLength - 1;
	my $downStart = $targetCandidate[2];
	my $downEnd = $targetCandidate[2] + $flankLength - 1;
	
	if ($upStart < 0)
	{
		$upStart = 0;
	}
	
	if ($downEnd >= length($rnaSeq))
	{
		$downEnd = length($rnaSeq) - 1;
	}
	
	my $upFlankSeq = substr($rnaSeq, $upStart,
	                        $upEnd - $upStart + 1);
	$upFlankSeq = reverse($upFlankSeq);
	my $downFlankSeq = substr($rnaSeq, $downStart, 
	                          $downEnd - $downStart + 1);
	                          
	my $upFlankAuCount = 0;
	my $downFlankAuCount = 0;
	
	# Need separate loops as flanks may have different lengths
	for (my $i = 0; $i < length($upFlankSeq); $i++)
	{
		my $letter = substr($upFlankSeq, $i, 1);
		
		if ($letter eq "A" || $letter eq "U")
		{
			$upFlankAuCount += 1 / ($i + 1);
		}
	}
	
	for (my $i = 0; $i < length($downFlankSeq); $i++)
	{
		my $letter = substr($downFlankSeq, $i, 1);
		
		if ($letter eq "A" || $letter eq "U")
		{
			$downFlankAuCount += 1 / ($i + 1);
		}
	}
	
	return($upFlankAuCount + $downFlankAuCount);
}

1;
