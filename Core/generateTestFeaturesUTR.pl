#!/usr/bin/perl

# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
# use warnings;

use Bio::SeqIO;
use List::Util qw(reduce);

require 'rna_comp_func.pl';
require 'rna_cons_func.pl';
require 'rna_match_func.pl';
require 'rna_misc_func.pl';
require 'rna_struct_func.pl';
require 'rna_miranda_func.pl';

my $mirFileName = shift();
my $utrFileName = shift();
my $phastFileName = shift();
my $pairsFileName = shift();
my $class = shift() || "?";
my $minTargetSiteLen = shift() || 10;
my $minMirandaScore = shift() || 80;
my %mirs = ();
my %utrs = ();
my %phastCons = ();

# Read in mature miRs from file
my $seqIO = Bio::SeqIO->new(-file => "<$mirFileName", 
                             -format => "fasta");

while (my $seq = $seqIO->next_seq) {
	$mirs{$seq->id} = $seq;
}

# Read in 3' UTRs from file
$seqIO = Bio::SeqIO->new(-file => "<$utrFileName", 
                          -format => "fasta");

while (my $seq = $seqIO->next_seq) {
	#my $id = $seq->id;
	#print STDERR "$id\n";
	#my $desc = $seq->description;
	#print STDERR "$desc\n";
	my @ids = ($seq->id =~ m/(N[MR]_\d+)/g);
	$utrs{$ids[0]} = $seq;
}

# Read in PhastConsScores
open(PHAST, $phastFileName) or die("Cannot open $phastFileName\n");

while (my $line = <PHAST>)
{
	chomp($line);
	my @tokens = split(/ /, $line);
	my @ids = ($tokens[0] =~ m/(N[MR]_\d+)_\.*/g);
	#print STDERR "$ids[0]\n";
	my @scores = @tokens[1 .. $#tokens];
	#print STDERR "@scores\n";
	$phastCons{$ids[0]} = \@scores;
}

close(PHAST);

# Print header
print("Miranda_score\t");
print("miR_ID\t");
print("mRNA_ID\t");
print("Start_position\t");
print("Duplex_MFE\t");
print("Seed_MFE\t");
print("3p_MFE\t");
print("miR_match_P04\t");
print("miR_match_P17\t");
print("miR_match_P20\t");
print("Seed_match_7mer1\t");
print("Seed_match_7mer2\t");
print("Seed_match_7merA1\t");
print("Seed_match_6mer2\t");
print("Seed_match_6mer1GU\t");
print("Seed_AU\t");
print("3p_mismatch\t");
print("Class\n");

# Process target site file
open(my $fh, "<", "$pairsFileName") or die("Cannot open $pairsFileName\n");
my @lines = readline($fh);
close($fh);

chomp(@lines);
my @listOfPairs = map { my @result = split(/\t+/, $_); \@result } @lines;
my %hash = ();
unshift(@listOfPairs, \%hash);
my $hashTable = reduce {
        my %localHash = %$a;
        my @localPair = @$b;
	    my @undefArray = ();
        my $ref = $localHash{$localPair[0]} || \@undefArray;
        my @currentArray = @$ref;
        push(@currentArray, $localPair[1]);
        $localHash{$localPair[0]} = \@currentArray;
	\%localHash
    } @listOfPairs;

my %utrHashTable = %$hashTable;

for my $utrId (keys %utrHashTable) {
    my @mirArray = @{$utrHashTable{$utrId}};

	if(!defined($utrs{$utrId})) {
        print STDERR "Missing utr: $utrId\n";
        next;
	}

	my $mrnaObj = $utrs{$utrId};
	my $mrnaSeq = formatRnaSeq($mrnaObj->seq);
	my @accMatrix = getAccessibility($mrnaSeq);
	my @phastConsScores = @{$phastCons{$utrId}};

	if (scalar(@phastConsScores) != length($mrnaObj->seq)) {
	    next;
	}

    for my $mirId (@mirArray) {

	    # Skip header
	    if ($mirId eq "miR_ID") {
	    	next;
	    }

	    if(!defined($mirs{$mirId})) {
            print STDERR "Missing mir: $mirId\n";
            next;
	    }

		print STDERR "$mirId $utrId...\n";
		my $mirObj = $mirs{$mirId};
		my $mirSeq = formatRnaSeq($mirObj->seq);

		my ($candidates, $mrnaAlign, $mirAlign) = findTargetCandidatesMiranda($mrnaObj->seq,  $mirObj->seq, $minMirandaScore);
		my @candidateTargets = @{$candidates};
		my @mrnaAlignments = @{$mrnaAlign};
		my @mirAlignments = @{$mirAlign};
		my $len = scalar(@candidateTargets);

		for (my $i = 0; $i < $len; $i++) {
			my @candidateTarget = @{$candidateTargets[$i]};
			my $startPos = $candidateTarget[1];
			my $endPos = $candidateTarget[2];
			my $mirandaScore = $candidateTarget[6];
			my @mrnaAlignment = @{$mrnaAlignments[$i]};
			my @mirAlignment = @{$mirAlignments[$i]};
			
			if ($endPos - $startPos + 1 < $minTargetSiteLen) {
				next;
			}
						
			while ($mirAlignment[-1] eq "-") {
				pop(@mrnaAlignment);
				pop(@mirAlignment);
			}
			
			if (scalar(@mrnaAlignment) == scalar(@mirAlignment)) {
				print("$mirandaScore\t");
				print("$mirId\t$utrId\t$startPos\t");
				
				@candidateTarget = getDuplexEnergy($mrnaSeq, $mirSeq, \@candidateTarget);

				my ($duplexMfe, $seedMfe, $threeMfe) =  getDuplexMfe($mrnaSeq, $mirSeq, \@candidateTarget);

				print("$duplexMfe\t$seedMfe\t$threeMfe\t");

				my ($mirPosMatch, $seedTypes, $seedRegCounts, $threeRegCounts, $totalRegCounts) = getTargetSiteMatchFeatures(\@mrnaAlignment, \@mirAlignment);
				my @mirPosMatchCast = @$mirPosMatch;
				my @seedTypesCast = @$seedTypes;
				my @seedRegCountsCast = @$seedRegCounts;
				my @threeRegCountsCast = @$threeRegCounts;
				my @totalRegCountsCast = @$totalRegCounts;

				print("$mirPosMatchCast[14]\t");
				print("$mirPosMatchCast[16]\t");
				print("$mirPosMatchCast[19]\t");

                print("$seedTypesCast[2]\t");
                print("$seedTypesCast[3]\t");
                print("$seedTypesCast[4]\t");
                print("$seedTypesCast[5]\t");
                print("$seedTypesCast[7]\t");
                print("$seedTypesCast[8]\t");

                print("$seedRegCountsCast[1]\t"); # Seed_AU

                print("$threeRegCountsCast[3]\t"); # 3p_mismatch

				print("$class\n");
			}
		}
	}
}
