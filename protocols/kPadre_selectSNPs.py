#!/usr/bin/env python3

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 
##                              kPadre selectSNPs 
##
##		This script parses a variant calls file containing variants
##		from parent 2 relative to genome of parent 1 and creates 
##		a SNPs file for "kPadre oneGenome".
##
##		input: variants call file in .vcf format 
##		output: SNPs file in .txt format 
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# import modules
import sys
import re
import subprocess

## define parsing function
def selectSNPs(vcfFile, bs):
	''' Parse variants in VCF file and selects diagnostic SNPs 
	
		-input: 
			vcfFile(str) = path to vcf.file
			ns (bolean)  = set True if processing bisulfite sequencing data 
		-output:
			SNPs file lines (stout)
	''' 

	with open(vcfFile, 'r') as fh_in:

		for line in fh_in.readlines():

			# strip formating
			line = line.strip().split()

			# Discard header lines
			if line[0].startswith('#'):
				continue
				
			# Obtain values
			chr = line[0]
			pos = line[1]
			ref = line[3]
			alts = [x for x in line[4].split(',')]
			desc = line[9].split(':')[0]
			geno = (desc[0], desc[2])

			# discard sites that contain reference allele
			# or no genotype

			if '0' in geno:
				continue
			if '.' in geno:
				continue
			if '*' in geno:
				continue
			
			# Obtain alternative allele 
			index = int(geno[0]) - 1
			alt_geno = alts[index]

			# Skip undetected indels
			if len(ref)>1 or len(alt_geno)>1:
				continue

			# Skip sites with undefined genotype
			if ref == '*' or alt_geno == '*':
				continue

			# Skip potential bisulphite transformation
			if bs:

				# Skip C -> T mutations for bisulfite transformation
				if ref == 'C' and alt_geno == 'T':
					continue

				# Skip G -> A mutations for bisulfite transformation in oposite strand
				if ref == 'G' and alt_geno == 'A':
					continue


			# Write line for position
			values = [chr, pos, ref, alt_geno]
			print(f"{' '.join(values)}")

# define main funciton
def main(vcfFile, bs):
	selectSNPs(vcfFile, bs)

# Function for dispatcher kPadre.py
def run(vcfFile, bs):
	main(vcfFile, bs)

# Run as stand alone
if __name__ == "__main__":

	# catch wrong syntax
	if len(sys.argv) > 3:
		print("Usage: ./kPadre_selectSNPs variants.vcf [bs]")
		sys.exit()

	# Obtain positional arguments
	vcfFile = sys.argv[1]
	print(sys.argv)
	bs = True if 'bs' in sys.argv else False

	main(vcfFile, bs)