# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

# Functions related to RNA secondary structure via Vienna RNA package 
# for use in miR target prediction
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
#use File::Temp qw(tempfile);

# Use RNAduplex to identify target candidates of mirSeq in mrnaSeq
sub findTargetCandidates($$$)
{
	my $mrnaSeq = shift();
	my $mirSeq = shift();
	my $eRange = shift() || 10;
	
	#my ($tempFH, $tempFN) = tempfile(UNLINK => 1);
	#my ($tempFH, $tempFN) = tempfile();
	#print("$tempFH\n$tempFN\n");
	#print($tempFH "$mrnaSeq\n$mirSeq\n");
	#$tempFH->flush();
	
	my @targetCandidates = ();
	#open(OUTPUT, "RNAduplex -s -e $eRange < $tempFN 2> /dev/null |");
	#print("$mrnaSeq\n");
	#print("echo -e \"$mrnaSeq\n$mirSeq\n\"");
	#print("echo \"$mrnaSeq\n$mirSeq\n\"");
	open(OUTPUT, "echo \"$mrnaSeq\n$mirSeq\n\" | RNAduplex -s -e $eRange 2> /dev/null |");
	#open(OUTPUT, "RNAduplex -s -e $eRange < test.txt 2> /dev/null |");
	#print(my $line = <OUTPUT>);
	my $i = 0;
	while (my $line = <OUTPUT>) 
	{
		chomp($line);
		#print("$line\n");
		#my @tokens = split(/\s+/, $line);
		#my @mrnaPosTokens = split(/,/, $tokens[1]);
		#my @mirPosTokens = split(/,/, $tokens[3]);
		#print("@tokens\n");
		#$targetCandidates[$i][0] = $tokens[0];
		#$targetCandidates[$i][1] = $mrnaPosTokens[0];
		#$targetCandidates[$i][2] = $mrnaPosTokens[1];
		#$targetCandidates[$i][3] = $mirPosTokens[0];
		#$targetCandidates[$i][4] = $mirPosTokens[1];
		#$targetCandidates[$i][5] = substr($tokens[4], 1, -1);
		
		$line =~ /(.*?)[ \s]+(\d+),(\d+)[ \s]+:[ \s]+(\d+),(\d+)[ \s]+\( *([-.\d]+)\)/;
		$targetCandidates[$i][0] = $1;
		$targetCandidates[$i][1] = $2;
		$targetCandidates[$i][2] = $3;
		$targetCandidates[$i][3] = $4;
		$targetCandidates[$i][4] = $5;
		$targetCandidates[$i][5] = $6;
		
		#for (my $j = 0; $j < 6; $j++) 
		#{
		#	print("$targetCandidates[$i][$j]\n");
		#}
		
		$i++;
	}
	
	close(OUTPUT);
	#close($tempFH);
	
	return(@targetCandidates);
}

# Convert bracket notation to an alignment
sub createAlignment($$\@)
{
	my $mrnaSeq = shift();
	my $mirSeq = shift();
	my @targetCandidate = @{shift()};
	
	#for (my $i = 0; $i < 6; $i++) 
	#{
	#	print("$targetCandidate[$i]\n");
	#}
	
	my @mrnaAlignment = ();
	my @mirAlignment = ();
	#print("$targetCandidate[1], $targetCandidate[2]\n");
	my $mrnaSiteSeq = substr($mrnaSeq, $targetCandidate[1]-1,
	                         $targetCandidate[2]-$targetCandidate[1]+1);
	my $numPrependSpaces = length($mirSeq)-$targetCandidate[4];
	my @mrnaSeqArray = split(//, $mrnaSiteSeq);
	my @mirSeqArray = split(//, $mirSeq);
	my @brackets = split(/&/, $targetCandidate[0]);
	my @mrnaBracket = split(//, $brackets[0]);
	my @mirBracket = split(//, $brackets[1]);
	#print("$mrnaSiteSeq\n");
	#print("$mirSeq\n");
	#print("$targetCandidate[0]\n");
	
	
	# Initial gaps due to extra miR sequence at 3' end
	for (my $i = 0; $i < $numPrependSpaces; $i++) 
	{
		push(@mrnaAlignment, "-");
		push(@mirAlignment, pop(@mirSeqArray));
	}
	
	my $left;
	while(defined($left = $mrnaBracket[0])) 
	{
		#print("$left\n");
		if ($left eq "." && $mirBracket[-1] eq ".") 
		{
			push(@mrnaAlignment, shift(@mrnaSeqArray));
			push(@mirAlignment, pop(@mirSeqArray));
			shift(@mrnaBracket);
			pop(@mirBracket);
			#print("@mrnaAlignment\n");
		}
		elsif ($left eq ".") 
		{
			push(@mrnaAlignment, shift(@mrnaSeqArray));
			push(@mirAlignment, "-");
			shift(@mrnaBracket);
		}
		elsif ($mirBracket[-1] eq ".") 
		{
			push(@mrnaAlignment, "-");
			push(@mirAlignment, pop(@mirSeqArray));
			pop(@mirBracket);
		}
		elsif ($left eq "(" && $mirBracket[-1] eq ")") 
		{
			push(@mrnaAlignment, shift(@mrnaSeqArray));
			push(@mirAlignment, pop(@mirSeqArray));
			shift(@mrnaBracket);
			pop(@mirBracket);
		}
		else 
		{
			# We got issues!
		}
	}
	
	# End gaps due to extra miR sequence at 5' end
	my $right;
	while(defined($right = pop(@mirSeqArray))) 
	{
		push(@mrnaAlignment, "-");
		push(@mirAlignment, $right);
	}
	
	
	#print("@mrnaAlignment\n");
	#print("@mirAlignment\n");
	
	#my @alignment;
	#@alignment[0] = @mrnaAlignment;
	#@alignment[1] = @mirAlignment;
	
	#print("@alignment[0]\n");
	#print("@alignment[1]\n");
	
	#return(@alignment);
	return(\@mrnaAlignment, \@mirAlignment);
}

# Get MFE of duplex of candidate target site
sub getDuplexMfe($$\@)
{
	my $mrnaSeq = shift();
	my $mirSeq = shift();
	my @targetCandidate = @{shift()};
	
	my $targetSite = substr($mrnaSeq, $targetCandidate[1]-1, 
	                        $targetCandidate[2]-$targetCandidate[1]+1);
	my $seedReg = substr($mirSeq, 0, 8);
	my $threeReg = substr($mirSeq, 8); #my $threeReg = substr($mirSeq, 8, -1);
	my $duplexMfe = $targetCandidate[5];
	my $seedMfe = 0;
	my $threeMfe = 0;
	
	open(OUTPUT, "echo \"$targetSite\n$seedReg\n\" | RNAduplex 2> /dev/null |");
	my $line = <OUTPUT>;
	chomp($line);		
	#print("$line\n");
	$line =~ /(.*?)[ \s]+(\d+),(\d+)[ \s]+:[ \s]+(\d+),(\d+)[ \s]+\( *([-.\d]+)\)/;
	$seedMfe = $6;
	
	open(OUTPUT, "echo \"$targetSite\n$threeReg\n\" | RNAduplex 2> /dev/null |");
	$line = <OUTPUT>;
	chomp($line);	
	#print("$line\n");	
	$line =~ /(.*?)[ \s]+(\d+),(\d+)[ \s]+:[ \s]+(\d+),(\d+)[ \s]+\( *([-.\d]+)\)/;
	$threeMfe = $6;
	
	return($duplexMfe, $seedMfe, $threeMfe);
}

# Return MFE of unconstrained RNA using RNAfold
sub getRnaMfe($)
{
	my $rnaSeq = shift();
	
	open(OUTPUT, "echo \"$rnaSeq\n\" | RNAfold --noPS 2> /dev/null |");

	#while (my $line = <OUTPUT>) 
	#{
	#	print("$line");
	#}
	
	# Get second line of output
	my $line = <OUTPUT>;
	$line = <OUTPUT>;
	close(OUTPUT);
	chomp($line);
	$line =~ /(.*?)[ \s]+\( *([-.\d]+)\)/;
	
	return($2);
}

# Return local MFE around a target site of unconstrained RNA
sub getLocalRnaMfe($\@$)
{
	my $rnaSeq = shift();
	my @targetCandidate = @{shift()};
	my $window = shift() || 100;
	
	my $startPos = $targetCandidate[1] - $window - 1;
	my $endPos = $targetCandidate[2] + $window - 1;
	
	if ($startPos < 0)
	{
		$startPos = 0;
	}
	
	if ($endPos >= length($rnaSeq))
	{
		$endPos = length($rnaSeq) - 1;
	}
	
	my $localRnaSeq = substr($rnaSeq, $startPos, $endPos-$startPos+1);
	
	open(OUTPUT, "echo \"$localRnaSeq\n\" | RNAfold --noPS 2> /dev/null |");

	# Get second line of output
	my $line = <OUTPUT>;
	$line = <OUTPUT>;
	close(OUTPUT);
	chomp($line);
	$line =~ /(.*?)[ \s]+\( *([-.\d]+)\)/;
	
	return($2);
}

# Return MFE of constrained RNA using RNAfold
sub getConstRnaMfe($\@)
{
	my $rnaSeq = shift();
	my @targetCandidate = @{shift()};
	
	my $startConst = "." x ($targetCandidate[1] - 1);
	my $targetConst = "x" x ($targetCandidate[2] - $targetCandidate[1] + 1);
	my $endConst = "." x (length($rnaSeq) - $targetCandidate[2]);
	my $const = $startConst . $targetConst . $endConst;
	
	open(OUTPUT, "echo \"$rnaSeq\n$const\n\" | RNAfold --noPS -C 2> /dev/null|");

	#while (my $line = <OUTPUT>) 
	#{
	#	print("$line");
	#}
	
	# Get second line of output
	my $line = <OUTPUT>;
	$line = <OUTPUT>;
	close(OUTPUT);
	chomp($line);
	$line =~ /(.*?)[ \s]+\( *([-.\d]+)\)/;
	
	return($2);
}

# Return local MFE around a target site of constrained RNA
sub getLocalConstRnaMfe($\@$)
{
	my $rnaSeq = shift();
	my @targetCandidate = @{shift()};
	my $window = shift() || 100;
	
	# Setup constraints
	my $startConst = "." x ($targetCandidate[1] - 1);
	my $targetConst = "x" x ($targetCandidate[2] - $targetCandidate[1] + 1);
	my $endConst = "." x (length($rnaSeq) - $targetCandidate[2]);
	my $const = $startConst . $targetConst . $endConst;
	
	# Setup local sequences
	my $startPos = $targetCandidate[1] - $window - 1;
	my $endPos = $targetCandidate[2] + $window - 1;
	
	if ($startPos < 0)
	{
		$startPos = 0;
	}
	
	if ($endPos >= length($rnaSeq))
	{
		$endPos = length($rnaSeq) - 1;
	}
	
	my $localRnaSeq = substr($rnaSeq, $startPos, $endPos-$startPos+1);
	my $localConst = substr($const, $startPos, $endPos-$startPos+1);
	
	open(OUTPUT, "echo \"$localRnaSeq\n$localConst\n\" | RNAfold --noPS -C 2> /dev/null|");
	
	# Get second line of output
	my $line = <OUTPUT>;
	$line = <OUTPUT>;
	close(OUTPUT);
	chomp($line);
	$line =~ /(.*?)[ \s]+\( *([-.\d]+)\)/;
	
	return($2);
}

# Return opening energy of site
sub getOpeningEnergy($$)
{
	my $unconstMfe = shift();
	my $constMfe = shift();
	
	return ($unconstMfe - $constMfe);
}

# Return accessibility matrix of RNA sequence
sub getAccessibility($$$$)
{	
	my $rnaSeq = shift();
	my $window = shift() || 80;
	my $span = shift() || 40;
	my $uLength = shift() || 10;
	
	my @accMatrix = ();
	
	# Call RNAplfold which should generate a file named plfold_lunp
	system("echo \"$rnaSeq\n\" | RNAplfold -W $window -L $span -u $uLength");
	open(OUTPUT, "plfold_lunp") or die("Cannot open plfold_lunp!");
	
	# Skip first header lines
	my $line = <OUTPUT>;
	$line = <OUTPUT>;
	
	# Read in matrix from file
	my $i = 0;
	while ($line = <OUTPUT>)
	{
		chomp($line);
		my @tokens = split(/\t/, $line);
		
		if ($tokens[0] != $i + 1)
		{
			die("Error reading plfold_lunp!");
		}
		
		for (my $j = 1; $j <= $uLength; $j++)
		{
			$accMatrix[$i][$j-1] = $tokens[$j];
		}
		
		$i++
	}
	close(OUTPUT);
	
	if ($i != length($rnaSeq))
	{
		die("Error reading plfold_lunp!");
	}
	
	return(@accMatrix);
}

# Return accessibility of seed and flanking regions
sub getSeedAndFlankAccessibility(\@\@$)
{
	my @targetCandidate = @{shift()};
	my @accMatrix = @{shift()};
	my $flankLength = shift() || 10; # TODO: look up a reasonable value to use in lit
	
	my $nrows = @accMatrix;
	my $ncols = @{$accMatrix[0]};
	my @seedAcc = ();
	my @seed1merAcc = ();
	my @seedFlankAcc = ();
	my @up1merAcc = ();
	my @down1merAcc = ();
	my $targetSiteEnd = $targetCandidate[2];
	
	# 8-mer and 4-mer accessibility of seed region
	# Roughly due to potential bulges
	$seedAcc[0] = $accMatrix[$targetSiteEnd-1][7];
	# Roughly the accessibility of the 5' end of seed (relative to miR)
	$seedAcc[1] = $accMatrix[$targetSiteEnd-1][3];
	# Roughly the accessibility of the 3' end of seed (relative to miR)
	$seedAcc[2] = $accMatrix[$targetSiteEnd-6][3];
	# Currently only looking at 10 nt flank
	$seedFlankAcc[0] = $accMatrix[$targetSiteEnd-9][9];
	$seedFlankAcc[1] = $accMatrix[$targetSiteEnd+9][9];
        if(!looks_like_number($seedFlankAcc[0])) {
		$seedFlankAcc[0] = 1;      
        }
	
        if(!looks_like_number($seedFlankAcc[1])) {
		$seedFlankAcc[1] = 1;      
        }

	# 1-mer accessibility of seed region
	for (my $i = 0; $i < 8; $i++)
	{
		$seed1merAcc[$i] = $accMatrix[$targetSiteEnd-$i-1][0];
	}
	
	# 1-mer accessibility of seed flanks
	for (my $i = 0; $i < $flankLength; $i++)
	{
		my $j = $targetSiteEnd - $i - 9;
		if ($j >= 0)
		{
			$up1merAcc[$i] = $accMatrix[$j][0];
		}
		else
		{
			$up1merAcc[$i] = 1;
		}
		
		$j = $targetSiteEnd + $i - 1;
		if ($j < $nrows)
		{
			$down1merAcc[$i] = $accMatrix[$j][0];
		}
		else
		{
			$down1merAcc[$i] = 1;
		}
	}
	
	return (\@seedAcc, \@seed1merAcc, \@seedFlankAcc, \@up1merAcc, 
	        \@down1merAcc);
}


1;
