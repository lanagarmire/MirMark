#!/usr/bin/perl

# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;

my $wekaFn = shift();
my $dataFn = shift();
my @mirIds;
my @mrnaIds;
my @starts;
my @ends;
my @mirandaScores;

# Read in site information
open(DATAIN, $dataFn) or die("Cannot open $dataFn\n");

my $count = 0;
while(my $line = <DATAIN>)
{
	chomp($line);
	my @tokens = split(/\t/, $line);
	#print STDERR "@tokens\n";
	$mirandaScores[$count] = $tokens[0];
	$mirIds[$count] = $tokens[1];
	$mrnaIds[$count] = $tokens[2];
	$starts[$count] = $tokens[3];
	$ends[$count] = $tokens[4];
	#print STDERR "$count\t$mirIds[$count]\t$mrnaIds[$count]\t$starts[$count]\t$ends[$count]\n";
	$count++;
}

close(DATAIN);

open(IN, $wekaFn) or die("Cannot open $wekaFn\n");
my $inTestSection = 0;

while (my $line = <IN>)
{
	chomp($line);
	
	if ($line =~ /predicted/)
	{
		# Start of test data results
		$inTestSection = 1;
	}
	elsif ($inTestSection == 1 && $line eq "")
	{
		# End of test data results
		last;
	}
	elsif ($inTestSection)
	{
		#print STDERR "$line\n";
		my @tokens = split(/\s+/, $line);
		#print STDERR "@tokens\n";
		my $id = $tokens[1];
		my $mirandaScore = $mirandaScores[$id];
		my $mirName = $mirIds[$id];
		my $mrnaName = $mrnaIds[$id];
		my $startPos = $starts[$id];
		my $endPos = $ends[$id];
		my @prob;
		if ($tokens[4] eq "+" || $tokens[4] eq "-")
		{
			@prob = ($tokens[5] =~ m/([\d\.]+)/g);
		}
		else
		{
			@prob = ($tokens[4] =~ m/([\d\.]+)/g);
		}
		#print STDERR "@prob\n";
		print("$mirName\t$mrnaName\t$startPos\t$endPos\t$mirandaScore\t$prob[0]\n")
	}
}
close(IN);

