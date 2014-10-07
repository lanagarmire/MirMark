#!/usr/bin/perl

# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;

my $fn = shift();

# Process target site file
open(IN, $fn) or die("Cannot open $fn\n");

while (my $line = <IN>)
{
	chomp($line);
	my @tokens = split(/\t/, $line);
	print(join(",",@tokens[3 .. $#tokens]) . "\n");
}

close(IN);
