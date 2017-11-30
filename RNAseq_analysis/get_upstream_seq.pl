#!/usr/bin/perl
## this script is used to extract upstream sequences given: the transcription start site (TSS), 
## promoter size (upstream and downstream of the TSS), the reference genome sequences
## note: currently, this script will only work for human data and expects reference sequence for the X, Y and MT chromosomes
##
## argument 1: a file containing the probe name, chomrosome, start, stop, strand, ID and foldchange (generated from filter_expressions_with_locations.pl)
## argument 2: the location of the genome files, 1 file per chromosome file must be named: *X.fa, where X is the chromosome #
## argument 3: the amount of upstream sequence in base pairs
## argument 4: the amount of downstream sequence in base pairs
##
## Kimberly MacKay Nov. 1, 2016
## contact: kimberly.mackay@usask.ca
##
## license: This work is licensed under the Creative Commons Attribution-NonCommercial-
## ShareAlike 3.0 Unported License. To view a copy of this license, visit 
## http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 
## PO Box 1866, Mountain View, CA 94042, US

use strict;
use warnings;

## check to ensure three argument were passed in
die "ERROR: must pass in four argumnets." if @ARGV != 4;

## get the input arguments
my $sequences_to_get = $ARGV[0];
my $directory_for_genome_files = $ARGV[1];
my $upstream = $ARGV[2];
my $downstream = $ARGV[2];

## open the input file
open INFILE, "$sequences_to_get" or die "ERROR: $sequences_to_get could not be opened";
chomp(my @ensemblIDs = <INFILE>);
close INFILE;

## define some variables
my ($chr1_seq, $chr2_seq, $chr3_seq, $chr4_seq, $chr5_seq, $chr6_seq, $chr7_seq, $chr8_seq, $chr9_seq, $chr10_seq, 
	$chr11_seq, $chr12_seq, $chr13_seq, $chr14_seq, $chr15_seq, $chr16_seq, $chr17_seq, $chr18_seq, $chr19_seq, $chr20_seq,
	$chr21_seq, $chr22_seq, $chrX_seq, $chrY_seq, $chrMT_seq, $promoter_seq) = ("");

## open and read in each chromosome file
open CHR1, "$directory_for_genome_files". "1.fa" or die "ERROR: 1.fa  could not be opened";
chomp(my @chr1 = <CHR1>);
close CHR1;

foreach my $seq_line (@chr1)
{
	if(!($seq_line =~ />/))
	{
		$chr1_seq = $chr1_seq . $seq_line;
	}
}

open CHR2, "$directory_for_genome_files". "2.fa" or die "ERROR: 2.fa  could not be opened";
chomp(my @chr2 = <CHR2>);
close CHR2;

foreach my $seq_line (@chr2)
{
	if(!($seq_line =~ />/))
	{
		$chr2_seq = $chr2_seq . $seq_line;
	}
}

open CHR3, "$directory_for_genome_files". "3.fa" or die "ERROR: 3.fa  could not be opened";
chomp(my @chr3 = <CHR3>);
close CHR3;

foreach my $seq_line (@chr3)
{
	if(!($seq_line =~ />/))
	{
		$chr3_seq = $chr3_seq . $seq_line;
	}
}

open CHR4, "$directory_for_genome_files". "4.fa" or die "ERROR: 4.fa  could not be opened";
chomp(my @chr4 = <CHR4>);
close CHR4;

foreach my $seq_line (@chr4)
{
	if(!($seq_line =~ />/))
	{
		$chr4_seq = $chr4_seq . $seq_line;
	}
}

open CHR5, "$directory_for_genome_files". "5.fa" or die "ERROR: 5.fa  could not be opened";
chomp(my @chr5 = <CHR5>);
close CHR5;

foreach my $seq_line (@chr5)
{
	if(!($seq_line =~ />/))
	{
		$chr5_seq = $chr5_seq . $seq_line;
	}
}

open CHR6, "$directory_for_genome_files". "6.fa" or die "ERROR: 6.fa  could not be opened";
chomp(my @chr6 = <CHR6>);
close CHR6;

foreach my $seq_line (@chr6)
{
	if(!($seq_line =~ />/))
	{
		$chr6_seq = $chr6_seq . $seq_line;
	}
}

open CHR7, "$directory_for_genome_files". "7.fa" or die "ERROR: 7.fa  could not be opened";
chomp(my @chr7 = <CHR7>);
close CHR7;

foreach my $seq_line (@chr7)
{
	if(!($seq_line =~ />/))
	{
		$chr7_seq = $chr7_seq . $seq_line;
	}
}

open CHR8, "$directory_for_genome_files". "8.fa" or die "ERROR: 8.fa  could not be opened";
chomp(my @chr8 = <CHR8>);
close CHR8;

foreach my $seq_line (@chr8)
{
	if(!($seq_line =~ />/))
	{
		$chr8_seq = $chr8_seq . $seq_line;
	}
}

open CHR9, "$directory_for_genome_files". "9.fa" or die "ERROR: 9.fa  could not be opened";
chomp(my @chr9 = <CHR9>);
close CHR9;

foreach my $seq_line (@chr9)
{
	if(!($seq_line =~ />/))
	{
		$chr9_seq = $chr9_seq . $seq_line;
	}
}

open CHR10, "$directory_for_genome_files". "10.fa" or die "ERROR: 10.fa  could not be opened";
chomp(my @chr10 = <CHR10>);
close CHR10;

foreach my $seq_line (@chr10)
{
	if(!($seq_line =~ />/))
	{
		$chr10_seq = $chr10_seq . $seq_line;
	}
}

open CHR11, "$directory_for_genome_files". "11.fa" or die "ERROR: 11.fa  could not be opened";
chomp(my @chr11 = <CHR11>);
close CHR11;

foreach my $seq_line (@chr11)
{
	if(!($seq_line =~ />/))
	{
		$chr11_seq = $chr11_seq . $seq_line;
	}
}

open CHR12, "$directory_for_genome_files". "12.fa" or die "ERROR: 12.fa  could not be opened";
chomp(my @chr12 = <CHR12>);
close CHR12;

foreach my $seq_line (@chr12)
{
	if(!($seq_line =~ />/))
	{
		$chr12_seq = $chr12_seq . $seq_line;
	}
}

open CHR13, "$directory_for_genome_files". "13.fa" or die "ERROR: 13.fa  could not be opened";
chomp(my @chr13 = <CHR13>);
close CHR13;

foreach my $seq_line (@chr13)
{
	if(!($seq_line =~ />/))
	{
		$chr13_seq = $chr13_seq . $seq_line;
	}
}

open CHR14, "$directory_for_genome_files". "14.fa" or die "ERROR: 14.fa  could not be opened";
chomp(my @chr14 = <CHR14>);
close CHR14;

foreach my $seq_line (@chr14)
{
	if(!($seq_line =~ />/))
	{
		$chr14_seq = $chr14_seq . $seq_line;
	}
}

open CHR15, "$directory_for_genome_files". "15.fa" or die "ERROR: 15.fa  could not be opened";
chomp(my @chr15 = <CHR15>);
close CHR15;

foreach my $seq_line (@chr15)
{
	if(!($seq_line =~ />/))
	{
		$chr15_seq = $chr15_seq . $seq_line;
	}
}

open CHR16, "$directory_for_genome_files". "16.fa" or die "ERROR: 16.fa  could not be opened";
chomp(my @chr16 = <CHR16>);
close CHR16;

foreach my $seq_line (@chr16)
{
	if(!($seq_line =~ />/))
	{
		$chr16_seq = $chr16_seq . $seq_line;
	}
}

open CHR17, "$directory_for_genome_files". "17.fa" or die "ERROR: 17.fa  could not be opened";
chomp(my @chr17 = <CHR17>);
close CHR17;

foreach my $seq_line (@chr17)
{
	if(!($seq_line =~ />/))
	{
		$chr17_seq = $chr17_seq . $seq_line;
	}
}

open CHR18, "$directory_for_genome_files". "18.fa" or die "ERROR: 18.fa  could not be opened";
chomp(my @chr18 = <CHR18>);
close CHR18;

foreach my $seq_line (@chr18)
{
	if(!($seq_line =~ />/))
	{
		$chr18_seq = $chr18_seq . $seq_line;
	}
}

open CHR19, "$directory_for_genome_files". "19.fa" or die "ERROR: 19.fa  could not be opened";
chomp(my @chr19 = <CHR19>);
close CHR19;

foreach my $seq_line (@chr19)
{
	if(!($seq_line =~ />/))
	{
		$chr19_seq = $chr19_seq . $seq_line;
	}
}

open CHR20, "$directory_for_genome_files". "20.fa" or die "ERROR: 20.fa  could not be opened";
chomp(my @chr20 = <CHR20>);
close CHR20;

foreach my $seq_line (@chr20)
{
	if(!($seq_line =~ />/))
	{
		$chr20_seq = $chr20_seq . $seq_line;
	}
}

open CHR21, "$directory_for_genome_files". "21.fa" or die "ERROR: 21.fa  could not be opened";
chomp(my @chr21 = <CHR21>);
close CHR21;

foreach my $seq_line (@chr21)
{
	if(!($seq_line =~ />/))
	{
		$chr21_seq = $chr21_seq . $seq_line;
	}
}

open CHR22, "$directory_for_genome_files". "22.fa" or die "ERROR: 22.fa  could not be opened";
chomp(my @chr22 = <CHR22>);
close CHR22;

foreach my $seq_line (@chr22)
{
	if(!($seq_line =~ />/))
	{
		$chr22_seq = $chr22_seq . $seq_line;
	}
}

open CHRX, "$directory_for_genome_files". "X.fa" or die "ERROR: X.fa  could not be opened";
chomp(my @chrX = <CHRX>);
close CHRX;

foreach my $seq_line (@chrX)
{
	if(!($seq_line =~ />/))
	{
		$chrX_seq = $chrX_seq . $seq_line;
	}
}

open CHRY, "$directory_for_genome_files". "Y.fa" or die "ERROR: Y.fa  could not be opened";
chomp(my @chrY = <CHRY>);
close CHRY;

foreach my $seq_line (@chrY)
{
	if(!($seq_line =~ />/))
	{
		$chrY_seq = $chrY_seq . $seq_line;
	}
}

open CHRMT, "$directory_for_genome_files". "MT.fa" or die "ERROR: MT.fa  could not be opened";
chomp(my @chrMT = <CHRMT>);
close CHRMT;

foreach my $seq_line (@chrMT)
{
	if(!($seq_line =~ />/))
	{
		$chrMT_seq = $chrMT_seq . $seq_line;
	}
}

## define some variables
my ($probe_name, $chromosome, $start, $stop, $strand, $id, $exp, $promoter_start, $length);
$length = chomp($upstream) + chomp($downstream); 

## for each line in the input file
for(my $i = 1; $i <= $#ensemblIDs; $i++)
{
	## clear the promoter_seq varaible
	$promoter_seq = "";
	
	## split the line into the relavent bits
	($probe_name, $chromosome, $start, $stop, $strand, $id, $exp) = split /\s+/, $ensemblIDs[$i];
	
	## identify the promoter start site
	$promoter_start = $start - $upstream;
	
	## extract the promoter sequence from the reference genome sequence
	$promoter_seq = substr($chr1_seq, $promoter_start, $length) if $chromosome =~ /^1$/;
	$promoter_seq = substr($chr2_seq, $promoter_start, $length) if $chromosome =~ /^2$/;
	$promoter_seq = substr($chr3_seq, $promoter_start, $length) if $chromosome =~ /^3$/;
	$promoter_seq = substr($chr4_seq, $promoter_start, $length) if $chromosome =~ /^4$/;
	$promoter_seq = substr($chr5_seq, $promoter_start, $length) if $chromosome =~ /^5$/;
	$promoter_seq = substr($chr6_seq, $promoter_start, $length) if $chromosome =~ /^6$/;
	$promoter_seq = substr($chr7_seq, $promoter_start, $length) if $chromosome =~ /^7$/;
	$promoter_seq = substr($chr8_seq, $promoter_start, $length) if $chromosome =~ /^8$/;
	$promoter_seq = substr($chr9_seq, $promoter_start, $length) if $chromosome =~ /^9$/;
	$promoter_seq = substr($chr10_seq, $promoter_start, $length) if $chromosome =~ /^10$/;
	$promoter_seq = substr($chr11_seq, $promoter_start, $length) if $chromosome =~ /^11$/;
	$promoter_seq = substr($chr12_seq, $promoter_start, $length) if $chromosome =~ /^12$/;
	$promoter_seq = substr($chr13_seq, $promoter_start, $length) if $chromosome =~ /^13$/;
	$promoter_seq = substr($chr14_seq, $promoter_start, $length) if $chromosome =~ /^14$/;
	$promoter_seq = substr($chr15_seq, $promoter_start, $length) if $chromosome =~ /^15$/;
	$promoter_seq = substr($chr16_seq, $promoter_start, $length) if $chromosome =~ /^16$/;
	$promoter_seq = substr($chr17_seq, $promoter_start, $length) if $chromosome =~ /^17$/;
	$promoter_seq = substr($chr18_seq, $promoter_start, $length) if $chromosome =~ /^18$/;
	$promoter_seq = substr($chr19_seq, $promoter_start, $length) if $chromosome =~ /^19$/;
	$promoter_seq = substr($chr20_seq, $promoter_start, $length) if $chromosome =~ /^20$/;
	$promoter_seq = substr($chr21_seq, $promoter_start, $length) if $chromosome =~ /^21$/;
	$promoter_seq = substr($chr22_seq, $promoter_start, $length) if $chromosome =~ /^22$/;
	$promoter_seq = substr($chrX_seq, $promoter_start, $length) if $chromosome =~ /^X$/;
	$promoter_seq = substr($chrY_seq, $promoter_start, $length) if $chromosome =~ /^Y$/;
	$promoter_seq = substr($chrMT_seq, $promoter_start, $length) if $chromosome =~ /^MT$/;
	
	## print the promoter sequence in fasta format
	print ">$id\n";
	print "$promoter_seq\n";
}
