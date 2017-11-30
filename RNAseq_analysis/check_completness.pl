#!/usr/bin/perl
## this script is used to check that a multi-fasta file contains all gene IDs from a given list
## note: gene IDs must be in the fasta headers for this script to work
##
## argument 1: a text file that contains a list of the expected Ensembl Gene IDs
## argument 2: a multi-fasta file that contains the sequences, with Ensembl Gene IDs as the fatsa header
##
## Kimberly MacKay June 4, 2014
## contact: kimberly.mackay@usask.ca
##
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, US

use strict;
use warnings;

## check to ensure two arguments are passed in
if(@ARGV != 2)
{
	die "ERROR: must pass in two arguments.";
}

## get the input arguments
my $id_list = $ARGV[0];
my $sequence_file = $ARGV[1];

## define some variables
my %promoter_seq;
my $complete = 1;

## open the text file with the gene IDs
open IDS, "$id_list" or die "ERROR: $id_list could not be opened";
chomp(my @ids = <IDS>);
close IDS;

## open the mutli-fasta file
open SEQ, "$sequence_file" or die "ERROR: $sequence_file could not be opened";
chomp(my @sequences = <SEQ>);
close SEQ;

## initialize the hash containing the sequences identified by their fasta headers
for(my $k = 0; $k <= $#sequences-1; $k +=2)
{
	$promoter_seq{$sequences[$k]} = $sequences[$k+1];
}

## look up the sequences for the ids of interest
foreach my $id (@ids)
{
	## if the gene IDs are present in the multi-fasta file
	if(!(exists $promoter_seq{">$id"}) || $promoter_seq{">$id"} =~ /not/)
	{
		## update the boolean flag
		$complete = 0;
	}
	## otherwise, they are missing from the multi-fasta file
	else
	{
		## print the missing ID
		print "$id\n";
		
		## update the boolean flag
		$complete = 1;
	}
}

## if all the gene IDs were not present in the multi-fasta file
if(!$complete)
{
	## print a warning
	print "WARNING: $id_list is incomplete\n";
}