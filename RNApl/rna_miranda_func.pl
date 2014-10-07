# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

# Functions related to RNA secondary structure via Vienna RNA package 
# for use in miR target prediction
use strict;
use warnings;
use File::Temp qw(tempfile);

require 'rna_misc_func.pl';

# Use miRanda to identify target candidates of mirSeq in mrnaSeq
sub findTargetCandidatesMiranda($$$$)
{
	my $mrnaSeq = shift();
	my $mirSeq = shift();
	my $minScore = shift() || 100;
	my @targetCandidates = ();
	my @mrnaAlignments = ();
	my @mirAlignments = ();
	
	# Write out sequences to temporary files
	my ($tempUtrFh, $tempUtrFn) = tempfile(UNLINK => 1);
	#print("$tempUtrFh\n$tempUtrFn\n");
	print($tempUtrFh ">UTR\n$mrnaSeq\n");
	$tempUtrFh->flush();
	my ($tempMirFh, $tempMirFn) = tempfile(UNLINK => 1);
	#print("$tempMirFh\n$tempMirFn\n");
	print($tempMirFh ">miR\n$mirSeq\n");
	$tempMirFh->flush();
	
	open(OUTPUT, "miranda $tempMirFn $tempUtrFn -noenergy -sc $minScore |");
	
	my $i = 0;
	while (my $line = <OUTPUT>) 
	{
		chomp($line);
		#print("$line\n");
		
		if (substr($line, 0, 4) eq ">miR")
		{
			my @tokens = split(/\t/, $line);
			#print("@tokens\n");
			
			#foreach my $token (@tokens)
			#{
			#	print("$token\n");
			#}
			
			my @mrnaPos = split(/\s+/, $tokens[5]);
			my @mirPos = split(/\s+/, $tokens[4]);
			my @targetCandidate = ();
			$targetCandidate[0] = ""; # No bracket notation
			$targetCandidate[1] = $mrnaPos[0]; # mRNA start
			$targetCandidate[2] = $mrnaPos[1]; # mRNA end
			$targetCandidate[3] = $mirPos[0]; # miR start
			$targetCandidate[4] = $mirPos[1]; # miR end
			$targetCandidate[5] = ""; # No energy
			$targetCandidate[6] = $tokens[2]; # miRanda score	
			push(@targetCandidates, \@targetCandidate);
		}
		elsif ($line =~ /Query:/)
		{
			my @tokens = split(/\s+/, $line);
			#print("@tokens\n");
			my $seq = $tokens[3];
			$seq = formatRnaSeq($seq);
			#print("$seq\n");
			my @mirAlignment = split(//, $seq);
			push(@mirAlignments, \@mirAlignment);
			#print("@mirAlignment\n");
		}
		elsif ($line =~ /Ref:/)
		{
			my @tokens = split(/\s+/, $line);
			#print("@tokens\n");
			my $seq = $tokens[3];
			$seq = formatRnaSeq($seq);
			#print("$seq\n");
			my @mrnaAlignment = split(//, $seq);
			push(@mrnaAlignments, \@mrnaAlignment);
			#print("@mrnaAlignment\n");
		}
	}
	
	close(OUTPUT);
	close($tempUtrFh);
	close($tempMirFh);
	
	return(\@targetCandidates, \@mrnaAlignments, \@mirAlignments);
}

# Use RNAduplex to get duplex energy
sub getDuplexEnergy($$\@)
{
	my $mrnaSeq = shift();
	my $mirSeq = shift();
	my @targetCandidate = @{shift()};
	my $mrnaStart = $targetCandidate[1];
	my $mrnaEnd = $targetCandidate[2];
	my $mirStart = $targetCandidate[3];
	my $mirEnd = $targetCandidate[4];
	$mrnaSeq = substr($mrnaSeq, $mrnaStart - 1, 
	                    $mrnaEnd - $mrnaStart + 1);
	$mirSeq = substr($mirSeq, $mirStart - 1, 
	                   $mirEnd - $mirStart + 1);
	
	open(OUTPUT, "echo \"$mrnaSeq\n$mirSeq\n\" | RNAduplex -s 2> /dev/null |");
	
	my $line = <OUTPUT>;
	if (defined($line)) 
	{
		chomp($line);
		#print("$line\n");		
		$line =~ /(.*?)[ \s]+(\d+),(\d+)[ \s]+:[ \s]+(\d+),(\d+)[ \s]+\( *([-.\d]+)\)/;
		$targetCandidate[5] = $6;
	}
	
	close(OUTPUT);
	
	return(@targetCandidate);
}

1;
