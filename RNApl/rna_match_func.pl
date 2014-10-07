# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;

# Get numerical value of the match status between two alignment 
# positions
sub getMatchStatus($$)
{
	my $t = shift();
	my $b = shift();
	
	my $status = 0;
	
	if (($t eq "G" && $b eq "C") || ($t eq "C" && $b eq "G")) 
	{
		$status = 1;
	}
	elsif (($t eq "A" && $b eq "U") || ($t eq "U" && $b eq "A")) 
	{
		$status = 2;
	}
	elsif (($t eq "G" && $b eq "U") || ($t eq "U" && $b eq "G")) 
	{
		$status = 3;
	}
	elsif ($t eq "-" || $b eq "-") 
	{
		$status = 5;
	}
	else 
	{
		$status = 4;
	}
	
	return($status);
}

# Return several array references to features relating to target site
# matching
sub getTargetSiteMatchFeatures(\@\@)
{
	my @mrnaAlignment = @{shift()};
	my @mirAlignment = @{shift()};
	
	#print("@mrnaAlignment\n");
	#print("@mirAlignment\n");
	
	# miR position match status
	my @mirPosMatch = ();
	
	# Seed match types
	my $sm8mer = 1; #z1-z8 WC match
	my $sm8merA1 = 1; #z1 match/mismatch to A, z2-z8 WC match
	my $sm7mer1 = 1; #z1-z7 WC match
	my $sm7mer2 = 1; #z2-z8 WC match
	my $sm7merA1 = 1; #z1 match/mismatch to A, z2-z7 WC match
	my $sm6mer1 = 1; # z1-z6 WC match
	my $sm6mer2 = 1; # z2-z7 WC match
	my $sm6mer1GU = 1; # z1-z6 WC match or GU, only one GU match
	my $sm6mer2GU = 1; # z2-z7 WC match or GU, only one GU match
	
	# Seed region counts
	#my $rgsMatch = 0;
	my $rgsGC = 0;
	my $rgsAU = 0;
	my $rgsGU = 0;
	my $rgsMismatch = 0;
	my $rgsBulge = 0;
	my $rgsBulgeNt = 0;
	
	# 3' region counts
	#my $rg3Match = 0;
	my $rg3GC = 0;
	my $rg3AU = 0;
	my $rg3GU = 0;
	my $rg3Mismatch = 0;
	my $rg3Bulge = 0;
	my $rg3BulgeNt = 0;
	
	# Total region counts
	#my $rgtMatch = 0;
	my $rgtGC = 0;
	my $rgtAU = 0;
	my $rgtGU = 0;
	my $rgtMismatch = 0;
	my $rgtBulge = 0;
	my $rgtBulgeNt = 0;
	
	# Compute miR match status and region counts
	my $count = 0;
	my $i = 1;
	
	while (defined($mirAlignment[-$i])) 
	{
		if ($mirAlignment[-$i] ne "-") 
		{
			# Get miR match status
			$mirPosMatch[$count] = 
			  getMatchStatus($mrnaAlignment[-$i], $mirAlignment[-$i]);
			  
			# WC match
			#if ($mirPosMatch[$count] == 1 || 
			#    $mirPosMatch[$count] == 2) 
			#{
			#	$rgtMatch++;
			#	if ($count >= 0 && $count < 8) 
			#	{
			#		$rgsMatch++;
			#	}
			#	else 
			#	{
			#		$rg3Match++;
			#	}
			#}
			
			# GC match
			if ($mirPosMatch[$count] == 1) 
			{
				$rgtGC++;
				if ($count >= 0 && $count < 8) 
				{
					$rgsGC++;
				}
				else 
				{
					$rg3GC++;
				}
			}
			
			# AU match
			if ($mirPosMatch[$count] == 2) 
			{
				$rgtAU++;
				if ($count >= 0 && $count < 8) 
				{
					$rgsAU++;
				}
				else 
				{
					$rg3AU++;
				}
			}
			
			# GU wobble
			if ($mirPosMatch[$count] == 3) 
			{
				$rgtGU++;
				if ($count >= 0 && $count < 8) 
				{
					$rgsGU++;
				}
				else 
				{
					$rg3GU++;
				}
			}
			
			# Mismatch
			if ($mirPosMatch[$count] == 4) 
			{
				$rgtMismatch++;
				if ($count >= 0 && $count < 8) 
				{
					$rgsMismatch++;
				}
				else 
				{
					$rg3Mismatch++;
				}
			}
			
			# Bulge on miR side
			if ($mirPosMatch[$count] == 5) 
			{
				$rgtBulgeNt++;
				if ($count >= 0 && $count < 8) 
				{
					$rgsBulgeNt++;
				}
				else 
				{
					$rg3BulgeNt++;
				}
				
				# Opening a new bulge
				if ($count == 0 || $mirPosMatch[$count-1] != 5) 
				{
					$rgtBulge++;
					if ($count >= 0 && $count < 8) 
					{
						$rgsBulge++;
					}
					else 
					{
						$rg3Bulge++;
					}
				}
			}
			
			$count++;
		}
		else 
		{
			# Bulge on mRNA side
			$rgtBulgeNt++;
			if ($count >= 0 && $count < 8) 
			{
				$rgsBulgeNt++;
			}
			else 
			{
				$rg3BulgeNt++;
			}
			
			# Opening a new bulge
			if ($i == 0 || $mirAlignment[-$i+1] ne "-") 
			{
				$rgtBulge++;
				if ($count >= 0 && $count < 8) 
				{
					$rgsBulge++;
				}
				else 
				{
					$rg3Bulge++;
				}
			}
		}
		
		# Check if first miR position doesn't match/mismatch to an A
		if ($count == 1 && $mrnaAlignment[-$i] ne "A") 
		{
			$sm8merA1 = 0;
			$sm7merA1 = 0;
		}
		
		$i++
	}
	
	# If miR is shorter than 20 nt, fill remaining features with 0s
	for (; $count < 20; $count++) 
	{
		$mirPosMatch[$count] = 0;
	}
		
	#print("@mirPosMatch\n");
	
	# sm8mer z1-z8 WC match
	# sm8merA1 z1 match/mismatch to A, z2-z8 WC match
	# sm7mer1 z1-z7 WC match
	# sm7mer2 z2-z8 WC match
	# sm7merA1 z1 match/mismatch to A, z2-z7 WC match
	# sm6mer1 z1-z6 WC match
	# sm6mer2 z2-z7 WC match
	# sm6mer1GU z1-z6 WC match or GU, only one GU match
	# sm6mer2GU z2-z7 WC match or GU, only one GU match
	
	# Compute seed type
	for ($i = 0; $i < 8; $i++) 
	{
		# Not WC match
		if ($mirPosMatch[$i] != 1 && $mirPosMatch[$i] != 2) 
		{
			# Mismatch in z1-z8
			if ($i >= 0 && $i < 8) 
			{
				$sm8mer = 0;
			}
			
			# Mismatch in z2-z8
			if ($i >= 1 && $i < 8)
			{
				$sm7mer2 = 0;
				$sm8merA1 = 0;
			}
			
			# Mismatch in z1-z7
			if ($i >= 0 && $i < 7) 
			{
				$sm7mer1 = 0;
			}
			
			# Mismatch in z1-z6
			if ($i >= 0 && $i < 6) 
			{
				$sm6mer1 = 0;
			}
			
			# Mismatch in z2-z7
			if ($i >= 1 && $i < 7) 
			{
				$sm6mer2 = 0;
				$sm7merA1 = 0;
			}
		}
	}
	
	# Check if we need to unset sm6mer1GU and sm6mer2GU
	if ($rgsGU > 3) 
	{
		# Too many GU wobbles in seed region
		$sm6mer1GU = 0;
		$sm6mer2GU = 0;
	}
	else 
	{
		my $nGU1 = 0;
		my $nGU2 = 0;
		for ($i = 0; $i < 8; $i++) 
		{
			# Not WC match and not GC wobble
			if ($mirPosMatch[$i] > 3) 
			{
				if ($i >= 0 && $i < 6) 
				{
					$sm6mer1GU = 0;
				}
				
				if ($i >= 1 && $i < 7) 
				{
					$sm6mer2GU = 0;
				}
			}
			# GC wobble
			elsif ($mirPosMatch[$i] == 3) 
			{
				if ($i >= 0 && $i < 6) 
				{
					$nGU1++;
				}
				
				if ($i >= 1 && $i < 7) 
				{
					$nGU2++;
				}
			}
		}
		
		# Too many GU wobbles in 6mer
		if ($nGU1 > 1) 
		{
			$sm6mer1GU = 0;
		}
		
		if ($nGU2 > 1) 
		{
			$sm6mer2GU = 0;
		}
	}
	
	my @seedTypes = ($sm8mer, $sm8merA1, $sm7mer1, $sm7mer2, $sm7merA1,
	                 $sm6mer1, $sm6mer2, $sm6mer1GU, $sm6mer2GU);
	my @seedRegCounts = ($rgsGC, $rgsAU, $rgsGU, $rgsMismatch, 
	                     $rgsBulge, $rgsBulgeNt);
	my @threeRegCounts = ($rg3GC, $rg3AU, $rg3GU, $rg3Mismatch, 
	                      $rg3Bulge, $rg3BulgeNt);
	my @totalRegCounts = ($rgtGC, $rgtAU, $rgtGU, $rgtMismatch, 
	                      $rgtBulge, $rgtBulgeNt);
	
	@mirPosMatch = @mirPosMatch[0 .. 19];
	
	return (\@mirPosMatch, \@seedTypes, \@seedRegCounts, 
	        \@threeRegCounts, \@totalRegCounts);
}

1;
