#!/usr/bin/perl

# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;

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
my $outputBaseName = shift();
my $modelFileName = shift();
my $class = shift() || "?";
my $minTargetSiteLen = shift() || 10;
my $minMirandaScore = shift() || 80;
my %mirs = ();
my %utrs = ();
my %phastCons = ();

# Read in mature miRs from file
my $seqIO = Bio::SeqIO->new(-file => "<$mirFileName",  -format => "fasta");
while (my $seq = $seqIO->next_seq) {
	$mirs{$seq->id} = $seq;
}

# Read in 3' UTRs from file
$seqIO = Bio::SeqIO->new(-file => "<$utrFileName",  -format => "fasta");
while (my $seq = $seqIO->next_seq) {
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
	my @scores = @tokens[1 .. $#tokens];
	$phastCons{$ids[0]} = \@scores;
}

close(PHAST);

my @headersFull = (
"Miranda_score",
"miR_ID",
"mRNA_ID",
"Start_position",
"End_position",
"Duplex_MFE",
"Seed_MFE",
"X3p_MFE",
"Local_target_MFE",
"Local_cons_target_MFE",
"Local_open_energy",
"miR_match_P01",
"miR_match_P02",
"miR_match_P03",
"miR_match_P04",
"miR_match_P05",
"miR_match_P06",
"miR_match_P07",
"miR_match_P08",
"miR_match_P09",
"miR_match_P10",
"miR_match_P11",
"miR_match_P12",
"miR_match_P13",
"miR_match_P14",
"miR_match_P15",
"miR_match_P16",
"miR_match_P17",
"miR_match_P18",
"miR_match_P19",
"miR_match_P20",
"Seed_match_8mer",
"Seed_match_8merA1",
"Seed_match_7mer1",
"Seed_match_7mer2",
"Seed_match_7merA1",
"Seed_match_6mer1",
"Seed_match_6mer2",
"Seed_match_6mer1GU",
"Seed_match_6mer2GU",
"Seed_GC",
"Seed_AU",
"Seed_GU",
"Seed_mismatch",
"Seed_bulge",
"Seed_bulge_nt",
"X3p_GC",
"X3p_AU",
"X3p_GU",
"X3p_mismatch",
"X3p_bulge",
"X3p_bulge_nut",
"Total_GC",
"Total_AU",
"Total_GU",
"Total_mismatch",
"Total_bulge",
"Total_bulge_nt",
"Seed_acc",
"Seed_5p_acc",
"Seed_P01_acc",
"Seed_P02_acc",
"Seed_P03_acc",
"Seed_P04_acc",
"Seed_P05_acc",
"Seed_P06_acc",
"Seed_P07_acc",
"Seed_P08_acc",
"Seed_3p_acc",
"Up_seed_flank_acc",
"Down_seed_flank_acc",
"Up_seed_P01_acc",
"Up_seed_P02_acc",
"Up_seed_P03_acc",
"Up_seed_P04_acc",
"Up_seed_P05_acc",
"Up_seed_P06_acc",
"Up_seed_P07_acc",
"Up_seed_P08_acc",
"Up_seed_P09_acc",
"Up_seed_P10",
"Down_seed_P01_acc",
"Down_seed_P02_acc",
"Down_seed_P03_acc",
"Down_seed_P04_acc",
"Down_seed_P05_acc",
"Down_seed_P06_acc",
"Down_seed_P07_acc",
"Down_seed_P08_acc",
"Down_seed_P09_acc",
"Down_seed_P10",
"Target_A_comp",
"Target_C_comp",
"Target_G_comp",
"Target_U_comp",
"Up_A_comp",
"Up_C_comp",
"Up_G_comp",
"Up_U_comp",
"Down_A_comp",
"Down_C_comp",
"Down_G_comp",
"Down_U_comp",
"Target_AA_comp",
"Target_AC_comp",
"Target_AG_comp",
"Target_AU_comp",
"Target_CA_comp",
"Target_CC_comp",
"Target_CG_comp",
"Target_CU_comp",
"Target_GA_comp",
"Target_GC_comp",
"Target_GG_comp",
"Target_GU_comp",
"Target_UA_comp",
"Target_UC_comp",
"Target_UG_comp",
"Target_UU_comp",
"Up_AA_comp",
"Up_AC_comp",
"Up_AG_comp",
"Up_AU_comp",
"Up_CA_comp",
"Up_CC_comp",
"Up_CG_comp",
"Up_CU_comp",
"Up_GA_comp",
"Up_GC_comp",
"Up_GG_comp",
"Up_GU_comp",
"Up_UA_comp",
"Up_UC_comp",
"Up_UG_comp",
"Up_UU_comp",
"Down_AA_comp",
"Down_AC_comp",
"Down_AG_comp",
"Down_AU_comp",
"Down_CA_comp",
"Down_CC_comp",
"Down_CG_comp",
"Down_CU_comp",
"Down_GA_comp",
"Down_GC_comp",
"Down_GG_comp",
"Down_GU_comp",
"Down_UA_comp",
"Down_UC_comp",
"Down_UG_comp",
"Down_UU_comp",
"Seed_flank_AU_score",
"Seed_cons_score",
"Target_cons_score",
"Flank_cons_score",
"Dist_to_end",
"classes"
);

my @headersTXT = (
"Miranda_score",
"miR_ID",
"mRNA_ID",
"Start_position",
"End_position",
"miR_match_P01",
"miR_match_P03",
"miR_match_P04",
"miR_match_P08",
"miR_match_P15",
"Seed_bulge",
"Total_AU",
"Total_mismatch",
"Total_bulge",
"Total_bulge_nt",
"Seed_P01_acc",
"Seed_cons_score"
);

my @headersCSV = (
"miR_match_P01",
"miR_match_P03",
"miR_match_P04",
"miR_match_P08",
"miR_match_P15",
"Seed_bulge",
"Total_AU",
"Total_mismatch",
"Total_bulge",
"Total_bulge_nt",
"Seed_P01_acc",
"Seed_cons_score"
);

open(TXT, ">$outputBaseName.txt");
open(CSV, ">$outputBaseName.csv");
open(ARFF, ">$outputBaseName.arff");

# Print headers on each file
print TXT join("\t", @headersTXT) . "\n";
print CSV "\"" . join("\",\"", @headersCSV) . "\"\n";
print ARFF '@relation ' . "$outputBaseName\n";
print ARFF "\n";
for my $csvHead (@headersCSV) {
    print ARFF '@attribute ' . $csvHead . " numeric\n";
}
print ARFF '@attribute classes {positive,negative}' . "\n";
print ARFF "\n";
print ARFF '@data' . "\n";

# Process target site file
open(my $fh, "<", "$pairsFileName") or die("Cannot open $pairsFileName\n");
my @lines = readline($fh);
close($fh);

chomp(@lines);
my @listOfPairs = map { my @result = split(/\t+/, $_); \@result } @lines;

my %utrHashTable = ();
for my $pair (@listOfPairs) {
    my @localPair = @$pair;
    my @currentArray = defined($utrHashTable{$localPair[0]}) ? @{$utrHashTable{$localPair[0]}} : ();
    push(@currentArray, $localPair[1]);
    $utrHashTable{$localPair[0]} = \@currentArray;
}

my @lookupHashes;

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
			    my %valHash = ();

                $valHash{'Miranda_score'}   = "$mirandaScore";
                $valHash{'miR_ID'}          = "$mirId";
                $valHash{'mRNA_ID'}         = "$utrId";
                $valHash{'Start_position'}  = "$startPos";
                $valHash{'End_position'}    = "$endPos";

				my ($mirPosMatch, $seedTypes, $seedRegCounts, $threeRegCounts, $totalRegCounts) = getTargetSiteMatchFeatures(\@mrnaAlignment, \@mirAlignment);
				my @mirPosMatchCast = @$mirPosMatch;
				my @seedTypesCast = @$seedTypes;
				my @seedRegCountsCast = @$seedRegCounts;
				my @threeRegCountsCast = @$threeRegCounts;
				my @totalRegCountsCast = @$totalRegCounts;

                $valHash{'miR_match_P01'}   = "$mirPosMatchCast[0]";
                $valHash{'miR_match_P03'}   = "$mirPosMatchCast[2]";
                $valHash{'miR_match_P04'}   = "$mirPosMatchCast[3]";
                $valHash{'miR_match_P08'}   = "$mirPosMatchCast[7]";
                $valHash{'miR_match_P15'}   = "$mirPosMatchCast[14]";

                $valHash{'Seed_bulge'}      = "$seedRegCountsCast[4]";

                $valHash{'Total_AU'}        = "$totalRegCountsCast[1]";
                $valHash{'Total_mismatch'}  = "$totalRegCountsCast[3]";
                $valHash{'Total_bulge'}     = "$totalRegCountsCast[4]";
                $valHash{'Total_bulge_nt'}  = "$totalRegCountsCast[5]";

				my ($seedAcc, $seed1merAcc, $seedFlankAcc, $up1merAcc, $down1merAcc) =  getSeedAndFlankAccessibility(\@candidateTarget, \@accMatrix);
				my @seed1merAccValCast = @$seed1merAcc;
                $valHash{'Seed_P01_acc'}    = "$seed1merAccValCast[0]";

				my ($seedScore, $siteScore, $flankScore) = getTargetSiteConservationFeatures(\@candidateTarget,  \@phastConsScores);
                $valHash{'Seed_cons_score'} = "$seedScore";

                print TXT join("\t", @valHash{@headersTXT}) . "\n";
                print CSV join(",", @valHash{@headersCSV}) . "\n";
                print ARFF join(",", @valHash{@headersCSV}) . ",?\n";

                push(@lookupHashes, \%valHash);
			}
		}
	}
}

close(TXT);
close(CSV);
close(ARFF);

# Execute classifier on ARFF
`java weka.classifiers.meta.FilteredClassifier -T $outputBaseName.arff -l $modelFileName -p 0 > $outputBaseName.weka`;

# Process target site file
open(my $weka, "<$outputBaseName.weka") or die("Cannot open $outputBaseName.weka\n");
my @wekaLines = readline($weka);
close($weka);
chomp(@wekaLines);

shift(@wekaLines);
shift(@wekaLines);
shift(@wekaLines);
shift(@wekaLines);
pop(@wekaLines);

my @wekaHeaders = split(/\s+/, shift(@wekaLines));

my @lookupHeaders = ( "miR_ID", "mRNA_ID", "Start_position", "End_position", "Miranda_score" );
my @resultHeaders = ( "miR_ID", "mRNA_ID", "Start_position", "End_position", "Miranda_score", "predicted", "probability" );
my @resultHeadersOutput = ( "MIR", "UTR", "Start position", "End position", "Miranda score", "Predicted", "Probability" );

open(RESULT, ">$outputBaseName.result");

print RESULT "\"" . join("\",\"", @resultHeadersOutput) . "\"\n";

map {
        my %hash;
        @hash{@wekaHeaders} = split(/\s+/, $_);
        my %lookupHash = %{$lookupHashes[$hash{'inst#'}-1]};
        @hash{@lookupHeaders} = @lookupHash{@lookupHeaders};
        if($hash{'predicted'} =~ /positive/) {
            $hash{'probability'} = $hash{'error'};
        } else {
            $hash{'probability'} = 1.0-$hash{'error'};
        }

        print RESULT join(",", @hash{@resultHeaders}) . "\n";
    } @wekaLines;

close(RESULT);

