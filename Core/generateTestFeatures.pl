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
print("Miranda_score\tmiR_ID\tmRNA_ID\tStart_position\t");
print("End_position\tDuplex_MFE\tSeed_MFE\t3p_MFE\t");
print("Local_target_MFE\tLocal_cons_target_MFE\tLocal_open_energy\t");
print("miR_match_P01\tmiR_match_P02\tmiR_match_P03\tmiR_match_P04\t");
print("miR_match_P05\tmiR_match_P06\tmiR_match_P07\tmiR_match_P08\t");
print("miR_match_P09\tmiR_match_P10\tmiR_match_P11\tmiR_match_P12\t");
print("miR_match_P13\tmiR_match_P14\tmiR_match_P15\tmiR_match_P16\t");
print("miR_match_P17\tmiR_match_P18\tmiR_match_P19\tmiR_match_P20\t");
print("Seed_match_8mer\tSeed_match_8merA1\tSeed_match_7mer1\t");
print("Seed_match_7mer2\tSeed_match_7merA1\tSeed_match_6mer1\t");
print("Seed_match_6mer2\tSeed_match_6mer1GU\tSeed_match_6mer2GU\t");
print("Seed_GC\tSeed_AU\tSeed_GU\tSeed_mismatch\tSeed_bulge\t");
print("Seed_bulge_nt\t3p_GC\t3p_AU\t3p_GU\t3p_mismatch\t3p_bulge\t");
print("3p_bulge_nut\tTotal_GC\tTotal_AU\tTotal_GU\tTotal_mismatch\t");
print("Total_bulge\tTotal_bulge_nt\tSeed_acc\tSeed_5p_acc\t");
print("Seed_P01_acc\tSeed_P02_acc\tSeed_P03_acc\tSeed_P04_acc\t");
print("Seed_P05_acc\tSeed_P06_acc\tSeed_P07_acc\tSeed_P08_acc\t");
print("Seed_3p_acc\tUp_seed_flank_acc\tDown_seed_flank_acc\t");
print("Up_seed_P01_acc\tUp_seed_P02_acc\tUp_seed_P03_acc\t");
print("Up_seed_P04_acc\tUp_seed_P05_acc\tUp_seed_P06_acc\t");
print("Up_seed_P07_acc\tUp_seed_P08_acc\tUp_seed_P09_acc\t");
print("Up_seed_P10\t");
print("Down_seed_P01_acc\tDown_seed_P02_acc\tDown_seed_P03_acc\t");
print("Down_seed_P04_acc\tDown_seed_P05_acc\tDown_seed_P06_acc\t");
print("Down_seed_P07_acc\tDown_seed_P08_acc\tDown_seed_P09_acc\t");
print("Down_seed_P10\t");
print("Target_A_comp\tTarget_C_comp\tTarget_G_comp\tTarget_U_comp\t");
print("Up_A_comp\tUp_C_comp\tUp_G_comp\tUp_U_comp\t");
print("Down_A_comp\tDown_C_comp\tDown_G_comp\tDown_U_comp\t");
print("Target_AA_comp\tTarget_AC_comp\tTarget_AG_comp\t");
print("Target_AU_comp\tTarget_CA_comp\tTarget_CC_comp\t");
print("Target_CG_comp\tTarget_CU_comp\tTarget_GA_comp\t");
print("Target_GC_comp\tTarget_GG_comp\tTarget_GU_comp\t");
print("Target_UA_comp\tTarget_UC_comp\tTarget_UG_comp\t");
print("Target_UU_comp\t");
print("Up_AA_comp\tUp_AC_comp\tUp_AG_comp\t");
print("Up_AU_comp\tUp_CA_comp\tUp_CC_comp\t");
print("Up_CG_comp\tUp_CU_comp\tUp_GA_comp\t");
print("Up_GC_comp\tUp_GG_comp\tUp_GU_comp\t");
print("Up_UA_comp\tUp_UC_comp\tUp_UG_comp\t");
print("Up_UU_comp\t");
print("Down_AA_comp\tDown_AC_comp\tDown_AG_comp\t");
print("Down_AU_comp\tDown_CA_comp\tDown_CC_comp\t");
print("Down_CG_comp\tDown_CU_comp\tDown_GA_comp\t");
print("Down_GC_comp\tDown_GG_comp\tDown_GU_comp\t");
print("Down_UA_comp\tDown_UC_comp\tDown_UG_comp\t");
print("Down_UU_comp\t");
print("Seed_flank_AU_score\tSeed_cons_score\tTarget_cons_score\t");
print("Flank_cons_score\tDist_to_end\tClass\n");

# Process target site file
open(IN, $pairsFileName) or die("Cannot open $pairsFileName\n");

while (my $line = <IN>)
{
	chomp($line);
	my @tokens = split(/\t/, $line);
	my $mirId = $tokens[1];
	my $utrId = $tokens[0];
	
	# Skip header
	if ($mirId eq "miR_ID")
	{
		next;
	}
	
	if (defined($mirs{$mirId}) && defined($utrs{$utrId}))
	{	
		print STDERR "$mirId $utrId...\n";
		my $mirObj = $mirs{$mirId};
		my $mrnaObj = $utrs{$utrId};
		my @phastConsScores = @{$phastCons{$utrId}};
		#print STDERR "@phastConsScores\n";
		
		if (scalar(@phastConsScores) != length($mrnaObj->seq))
		{
			next;
		}
		
		#print STDERR "$mirId $utrId\n";
		
		my ($candidates, $mrnaAlign, $mirAlign) = 
				findTargetCandidatesMiranda($mrnaObj->seq, 
											$mirObj->seq,
											$minMirandaScore);
		my @candidateTargets = @{$candidates};
		my @mrnaAlignments = @{$mrnaAlign};
		my @mirAlignments = @{$mirAlign};
		my $len = scalar(@candidateTargets);
		#print STDERR "$len\n";
		
		for (my $i = 0; $i < $len; $i++)
		{
			my @candidateTarget = @{$candidateTargets[$i]};
			my $startPos = $candidateTarget[1];
			my $endPos = $candidateTarget[2];
			my $mirandaScore = $candidateTarget[6];
			my @mrnaAlignment = @{$mrnaAlignments[$i]};
			my @mirAlignment = @{$mirAlignments[$i]};
			
			# Do I need these variables?
			my $targetLength = $endPos - $startPos + 1;
			my $desc = $mrnaObj->description;
			
			if ($targetLength < $minTargetSiteLen)
			{
				next;
			}
						
			while ($mirAlignment[-1] eq "-")
			{
				pop(@mrnaAlignment);
				pop(@mirAlignment);
			}
			
			if (scalar(@mrnaAlignment) == scalar(@mirAlignment))
			{
				print("$mirandaScore\t");
				print("$mirId\t$utrId\t$startPos\t$endPos\t");
				
				my $mirSeq = formatRnaSeq($mirObj->seq);
				my $mrnaSeq = formatRnaSeq($mrnaObj->seq);
				@candidateTarget = getDuplexEnergy($mrnaSeq, $mirSeq, 
												   \@candidateTarget);
				#print STDERR "@candidateTarget\n";
				
				my ($duplexMfe, $seedMfe, $threeMfe) = 
					getDuplexMfe($mrnaSeq, $mirSeq, \@candidateTarget);
				#print STDERR "$duplexMfe $seedMfe $threeMfe\n";
				
				print("$duplexMfe\t$seedMfe\t$threeMfe\t");
				
				my $mrnaLocalMfe = getLocalRnaMfe($mrnaSeq, 
													\@candidateTarget);
				#print STDERR "UTR Local MFE $mrnaLocalMfe\n";
				
				print("$mrnaLocalMfe\t");

				my $mrnaLocalConstMfe = getLocalConstRnaMfe($mrnaSeq, 
														\@candidateTarget);
				#print STDERR "UTR Const MFE $mrnaLocalConstMfe\n";
				
				print("$mrnaLocalConstMfe\t");
				
				my $localOpeningEnergy = 
					getOpeningEnergy($mrnaLocalMfe, $mrnaLocalConstMfe);
				#print STDERR "Local opening energy $localOpeningEnergy\n";
				
				print("$localOpeningEnergy\t");
				
				my ($mirPosMatch, $seedTypes, $seedRegCounts, 
					$threeRegCounts, $totalRegCounts) = 
					getTargetSiteMatchFeatures(\@mrnaAlignment, 
											   \@mirAlignment);
				#print STDERR "miR position match status: @$mirPosMatch\n";
				#print STDERR "Seed type: @$seedTypes\n";
				#print STDERR "Seed region counts: @$seedRegCounts\n";
				#print STDERR "3' region counts: @$threeRegCounts\n";
				#print STDERR "Total region counts: @$totalRegCounts\n";
				
				foreach my $match (@$mirPosMatch)
				{
					print("$match\t");
				}
				
				foreach my $type (@$seedTypes)
				{
					print("$type\t");
				}
				
				foreach my $seedRegCount (@$seedRegCounts)
				{
					print("$seedRegCount\t");
				}
				
				foreach my $threeRegCount (@$threeRegCounts)
				{
					print("$threeRegCount\t");
				}
				
				foreach my $totalRegCount (@$totalRegCounts)
				{
					print("$totalRegCount\t");
				}
				
				my @accMatrix = getAccessibility($mrnaSeq);
				my ($seedAcc, $seed1merAcc, $seedFlankAcc, $up1merAcc,
					$down1merAcc) = 
					getSeedAndFlankAccessibility(\@candidateTarget, 
												 \@accMatrix);
				#print STDERR "Seed Acc @$seedAcc\n";
				#print STDERR "Seed 1-mer Acc @$seed1merAcc\n";
				#print STDERR "Seed Flank Acc @$seedFlankAcc\n";
				#print STDERR "Upstream 1-mer Acc @$up1merAcc\n";
				#print STDERR "Downstream 1-mer Acc @$down1merAcc\n";
				
				foreach my $seedAccVal (@$seedAcc)
				{
					print("$seedAccVal\t");
				}
				
				foreach my $seed1merAccVal (@$seed1merAcc)
				{
					print("$seed1merAccVal\t");
				}
				
				foreach my $seedFlankAccVal (@$seedFlankAcc)
				{
					print("$seedFlankAccVal\t");
				}
				
				foreach my $up1merAccVal (@$up1merAcc)
				{
					print("$up1merAccVal\t");
				}
				
				foreach my $down1merAccVal (@$down1merAcc)
				{
					print("$down1merAccVal\t");
				}
				
				my ($targetMonomersHash, $upFlankMonomersHash, 
					$downFlankMonomersHash, $targetDimersHash, 
					$upFlankDimersHash, $downFlankDimersHash) = 
					getCompositionFeatures($mrnaSeq, \@candidateTarget);

				my $targetMonomers = 
					monomerHash2Array($targetMonomersHash);
				#print STDERR "Target monomers: @$targetMonomers\n";
				
				foreach my $targetMonomer (@$targetMonomers)
				{
					print("$targetMonomer\t");
				}

				my $upFlankMonomers = 
					monomerHash2Array($upFlankMonomersHash);
				#print STDERR "Upstream target flank monomers: @$upFlankMonomers\n";
				
				foreach my $upFlankMonomer (@$upFlankMonomers)
				{
					print("$upFlankMonomer\t");
				}

				my $downFlankMonomers = 
					monomerHash2Array($downFlankMonomersHash);
				#print STDERR "Downstream target flank monomers: @$downFlankMonomers\n";
				
				foreach my $downFlankMonomer (@$downFlankMonomers)
				{
					print("$downFlankMonomer\t");
				}

				my $targetDimers = dimerHash2Array($targetDimersHash);
				#print STDERR "Target dimers: @$targetDimers\n";
				
				foreach my $targetDimer (@$targetDimers)
				{
					print("$targetDimer\t");
				}

				my $upFlankDimers = dimerHash2Array($upFlankDimersHash);
				#print STDERR "Upstream target flank dimers: @$upFlankDimers\n";
				
				foreach my $upFlankDimer (@$upFlankDimers)
				{
					print("$upFlankDimer\t");
				}

				my $downFlankDimers = 
					dimerHash2Array($downFlankDimersHash);
				#print STDERR "Downstream target flank dimers: @$downFlankDimers\n";
				
				foreach my $downFlankDimer (@$downFlankDimers)
				{
					print("$downFlankDimer\t");
				}
				
				my $seedFlankAuContent = 
					getSeedFlankAUContent($mrnaSeq, \@candidateTarget);
				#print STDERR "Seed flank weighted AU content: $seedFlankAuContent\n";
				
				print("$seedFlankAuContent\t");

				my ($seedScore, $siteScore, $flankScore) = 
					getTargetSiteConservationFeatures(\@candidateTarget, 
													  \@phastConsScores);
				#print STDERR "$seedScore $siteScore $flankScore\n";
				
				print("$seedScore\t$siteScore\t$flankScore\t");
				
				my $distEnd = getTargetSiteEndDistance($mrnaSeq, 
					\@candidateTarget);
				print("$distEnd\t");
				
				print("$class\n");
			}
		}
	}
	else
	{
		print STDERR "Missing $mirId - $utrId\n";
	}
}

close(IN);
