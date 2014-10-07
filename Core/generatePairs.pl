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
my $pairsFileName = shift();
my %mirs = ();
my %utrs = ();

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

open(TXT, ">$pairsFileName");

foreach my $mirId (keys %mirs) {
	foreach my $utrId (keys %utrs) {
		print TXT "$utrId\t$mirId\n";
	}
}

close(TXT);
