#!/usr/bin/perl

# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;

use Bio::SeqIO;

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
print("miR_match_P01\t");
print("miR_match_P02\t");
print("miR_match_P03\t");
print("miR_match_P04\t");
print("miR_match_P05\t");
print("miR_match_P06\t");
print("miR_match_P07\t");
print("miR_match_P08\t");
print("miR_match_P09\t");
print("miR_match_P10\t");
print("miR_match_P11\t");
print("miR_match_P12\t");
print("miR_match_P13\t");
print("miR_match_P14\t");
print("miR_match_P15\t");
print("miR_match_P16\t");
print("miR_match_P17\t");
print("miR_match_P18\t");
print("miR_match_P19\t");
print("miR_match_P20\t");
print("Seed_match_7mer1\t");
print("Seed_match_7mer2\t");
print("Seed_match_7merA1\t");
print("Seed_match_6mer2\t");
print("Seed_match_6mer1GU\t");
print("Seed_AU\t");
print("Seed_bulge\t");
print("3p_mismatch\t");
print("Total_AU\t");
print("Total_mismatch\t");
print("Total_bulge\t");
print("Total_bulge_nt\t");
print("Seed_P01_acc\t");
print("Seed_cons_score\t");
print("Class\n");

# Process target site file
open(IN, $pairsFileName) or die("Cannot open $pairsFileName\n");

while (my $line = <IN>)
{
	chomp($line);
	my @tokens = split(/\t/, $line);
	my $mirId = $tokens[1];
	my $utrId = $tokens[0];
	
	# Skip header
	if ($mirId eq "miR_ID") {
		next;
	}
	
	if (defined($mirs{$mirId}) && defined($utrs{$utrId})) {
		print STDERR "$mirId $utrId...\n";
		my $mirObj = $mirs{$mirId};
		my $mrnaObj = $utrs{$utrId};
		my @phastConsScores = @{$phastCons{$utrId}};
		#print STDERR "@phastConsScores\n";
		
		if (scalar(@phastConsScores) != length($mrnaObj->seq)) {
			next;
		}
		
		my ($candidates, $mrnaAlign, $mirAlign) =
				findTargetCandidatesMiranda($mrnaObj->seq, 
											$mirObj->seq,
											$minMirandaScore);
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
				
				my $mirSeq = formatRnaSeq($mirObj->seq);
				my $mrnaSeq = formatRnaSeq($mrnaObj->seq);
				@candidateTarget = getDuplexEnergy($mrnaSeq, $mirSeq, 
												   \@candidateTarget);

				my ($duplexMfe, $seedMfe, $threeMfe) = 
					getDuplexMfe($mrnaSeq, $mirSeq, \@candidateTarget);

				print("$duplexMfe\t$seedMfe\t$threeMfe\t");

				my ($mirPosMatch, $seedTypes, $seedRegCounts, 
					$threeRegCounts, $totalRegCounts) = 
					getTargetSiteMatchFeatures(\@mrnaAlignment, 
											   \@mirAlignment);
				my @mirPosMatchCast = @$mirPosMatch;
				my @seedTypesCast = @$seedTypes;
				my @seedRegCountsCast = @$seedRegCounts;
				my @threeRegCountsCast = @$threeRegCounts;
				my @totalRegCountsCast = @$totalRegCounts;

				print("$mirPosMatchCast[0]\t");
				print("$mirPosMatchCast[2]\t");
				print("$mirPosMatchCast[3]\t");
				print("$mirPosMatchCast[7]\t");
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
                print("$seedRegCountsCast[4]\t"); # Seed_bulge

                print("$threeRegCountsCast[3]\t"); # 3p_mismatch

                print("$totalRegCountsCast[1]\t"); # total_au
                print("$totalRegCountsCast[3]\t"); # total_mismatch
                print("$totalRegCountsCast[4]\t"); # total_bulge
                print("$totalRegCountsCast[5]\t"); # total_bulge_nt

				my @accMatrix = getAccessibility($mrnaSeq);
				my ($seedAcc, $seed1merAcc, $seedFlankAcc, $up1merAcc, $down1merAcc) =  getSeedAndFlankAccessibility(\@candidateTarget, \@accMatrix);

				my @seed1merAccValCast = @$seed1merAcc;

                print("$seed1merAccValCast[0]\t"); # Seed_P01_acc

				my ($seedScore, $siteScore, $flankScore) =
					getTargetSiteConservationFeatures(\@candidateTarget, 
													  \@phastConsScores);

				print("$seedScore\t");
				
				print("$class\n");
			}
		}
	}
	else {
		print STDERR "Missing $mirId - $utrId\n";
	}
}

close(IN);
