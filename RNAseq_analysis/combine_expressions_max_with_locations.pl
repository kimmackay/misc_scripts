#!/usr/bin/perl
## this script takes a file containing the following columns form a seqMonk probe report:
##	- probe name
##	- chromosome
##	- start
##	- end
##	- strand
##	- ID
##	- treatment 1 foldchange
##	- treatment 2 foldchange
## and combines transcripts that belong to the same gene (as defined by the ID column)
## note: the input file MUST be sorted according to the ID column
## note: currently it only works for two treatments at a time
## note: this code is not elegant; it could be improved to reduce redundancies and make it more general
##
## argument 1: a text file containing the probe name, chromosome, start, end, strand, ID, treatment 1 foldchange, and treatment 2 foldchange sorted by gene id
## argument 2: the name of the output file
##
## Kimberly MacKay July 7, 2014
##
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, USA.

use strict;
use warnings;
use List::Util qw[min max];

## check to ensure two arguments were passed in
die "ERROR: must pass in two argumnets." if @ARGV != 2;

## get the input arguments
my $info_file = $ARGV[0];
my $output_file = $ARGV[1];

## open the input file
open INFO_FILE, "$info_file" or die "ERROR: $info_file could not be opened.";
chomp(my @info_file = <INFO_FILE>);
close INFO_FILE;

## open the output file
open OUTFILE, ">$output_file" or die "ERROR: $output_file could not be opened.";
printf OUTFILE "%-15s %10s %15s %15s %10s %20s %20s %20s", "Probe_Name", "Chromosome", "Start", "Stop", "Strand", "ID", "Treatment_1_Fold-Change", "Treatment_2_Fold-Change\n";

## define some variables
my ($probe_name, $chromosome, $start, $stop, $strand, $id,$T1_exp, $T2_exp);
my ($T1_exp_max, $T2_exp_max, $T1_exp_min, $T2_exp_min, $start_min);
my ($probe_name_next, $chromosome_next, $start_next, $stop_next, $strand_next, $id_next, $T1_exp_next, $T2_exp_next, $T2_exp_best, $T1_exp_best);
my $id_prev = "";

my @T1_expressions;
my @T2_expressions;

my @start_sites;
my $i;

## for each line in the input file
for($i =1; $i <= $#info_file; $i++)
{
	## split the line into the relavent bits
	($probe_name, $chromosome, $start, $stop, $strand, $id, $T1_exp, $T2_exp) = split /\s+/, $info_file[$i];
	
	## add the fold-changes and start sites into the corresponding collections (will be used for determining max/min later)
	push(@T1_expressions, $T1_exp);
	push(@T2_expressions, $T2_exp);
	push(@start_sites, $start);
	
	## if it is not the last line
	if($i != $#info_file)
	{
		## take a look at the information on the next line
		($probe_name_next, $chromosome_next, $start_next, $stop_next, $strand_next, $id_next, $T1_exp_next, $T2_exp_next) = split /\s+/, $info_file[$i+1];
	
		## if the information on the next line corresponds to a different gene
		if($id !~ /$id_next/)
		{
			## get the max and min fold-change for each treatment and the min start site
			$T2_exp_max = max(@T2_expressions);
			$T1_exp_max = max(@T1_expressions);
			
			$T2_exp_min = min(@T2_expressions);
			$T1_exp_min = min(@T1_expressions);
			
			$start_min = min(@start_sites);
			
			## determine which fold-change (max or min) corresponds to the largest fold-change -> this is what will be reported
			if($T2_exp_max > 1/$T2_exp_min)
			{
				$T2_exp_best = $T2_exp_max;
			}
			else
			{
				$T2_exp_best = $T2_exp_min;
			}
			
			if($T1_exp_max > 1/$T1_exp_min)
			{
				$T1_exp_best = $T1_exp_max;
			}
			else
			{
				$T1_exp_best = $T1_exp_min;
			}
			
			## print the information based for the gene with the largest fold-change for each treatment and the min start site
			printf OUTFILE "%-15s %10s %15s %15s %10s %20s %20s %20s", "$probe_name", "$chromosome", "$start_min", "$stop", "$strand", "$id", "$T1_exp_best", "$T2_exp_best\n";
			
			## clear the collections in preparation for the next gene
			undef(@T2_expressions);
			undef(@T1_expressions);
			undef(@start_sites);	
		}
	}
	## otherwise; it is the last line of the file and therefore the last gene in the file
	else 
	{
		## get the max and min fold-change for each treatment and the min start site	
		$T2_exp_max = max(@T2_expressions);
		$T1_exp_max = max(@T1_expressions);
		
		$T2_exp_min = min(@T2_expressions);
		$T1_exp_min = min(@T1_expressions);
		
		$start_min = min(@start_sites);
		
		## determine which fold-change (max or min) corresponds to the largest fold-change -> this is what will be reported
		if($T2_exp_max > 1/$T2_exp_min)
		{
			$T2_exp_best = $T2_exp_max;
		}
		else
		{
			$T2_exp_best = $T2_exp_min;
		}
		
		if($T1_exp_max > 1/$T1_exp_min)
		{
			$T1_exp_best = $T1_exp_max;
		}
		else
		{
			$T1_exp_best = $T1_exp_min;
		}
		
		## print the information based for the gene with the largest fold-change for each treatment and the min start site
		printf OUTFILE "%-15s %10s %15s %15s %10s %20s %20s %20s", "$probe_name", "$chromosome", "$start_min", "$stop", "$strand",  "$id", "$T1_exp_best", "$T2_exp_best\n";
	}
}