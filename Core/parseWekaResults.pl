#!/usr/bin/perl

# Copyright 2014 Lana Garmire, lgarmire@cc.hawaii.edu

use strict;
use warnings;

my $resDir = shift();

opendir(DIR, $resDir) or die $!;

print("NumFeatures,TPRate_P,FPRate_P,Prec_P,Rec_P,F_P,");
print("AUC_P,TPRate_N,FPRate_N,Prec_N,Rec_N,F_N,AUC_N,");
print("TP,FN,FP,TN\n");

while (my $file = readdir(DIR)) 
{
	next if ($file =~ m/^\./);
	next if ($file =~ m/.*\.model/);
	
	if ($file =~ m/10fold/)
	{
		#print STDERR "$file\n";
		
		open(IN, $resDir."/".$file) or die $!;
		my $accFlag = 0;
		my $conFlag = 0;
		my $cvFlag = 0;
		my $count = 0;
		my @numFeatures = ($file =~ m/top(\d+)\./);
		#print STDERR "@numFeatures\n";

		print("$numFeatures[0],");
		
		while (my $line = <IN>)
		{
			chomp($line);
			
			if ($conFlag == 1)
			{
				$count++;
				
				if ($count == 3)
				{
					#print STDERR "$line\n";
					my @tokens = split(/\s+/, $line);
					#print STDERR "@tokens";
					print("$tokens[1],$tokens[2],");
				}
				elsif ($count == 4)
				{
					#print STDERR "$line\n";
					my @tokens = split(/\s+/, $line);
					#print STDERR "@tokens";
					print("$tokens[1],$tokens[2]");
					print("\n");
				}
			}
			elsif ($accFlag == 1)
			{
				$count++;
				
				if ($count == 3)
				{
					#print STDERR "$line\n";
					my @tokens = split(/\s+/, $line);
					#print STDERR "@tokens";
					print("$tokens[1],$tokens[2],$tokens[3],");
					print("$tokens[4],$tokens[5],$tokens[6],");
				}
				elsif ($count == 4)
				{
					#print STDERR "$line\n";
					my @tokens = split(/\s+/, $line);
					#print STDERR "@tokens";
					print("$tokens[1],$tokens[2],$tokens[3],");
					print("$tokens[4],$tokens[5],$tokens[6],");
				}
			}
			
			if ($line eq "=== Stratified cross-validation ===")
			{
				$cvFlag = 1;
				next;
			}
			
			if ($cvFlag == 1 && $line eq "=== Detailed Accuracy By Class ===")
			{
				$accFlag = 1;
				next;
			}
			
			if ($accFlag == 1 && $line eq "=== Confusion Matrix ===")
			{
				$conFlag = 1;
				$count = 0;
				next;
			}
			
		}
		
		close(IN);
	}
}

closedir(DIR);
