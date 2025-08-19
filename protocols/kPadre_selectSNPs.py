## -------------------------------------------------------
## QueSNP: Scan vcf of filtered genotypes to create
##	   a table of genotypes of parental species
## -------------------------------------------------------

## -------------------------------------------------------
## Set up environment
## -------------------------------------------------------

# import modules .........................................
import sys
import re
import subprocess

## -------------------------------------------------------
## Obtain input and output files
## -------------------------------------------------------
vcf = sys.argv[1]
out = sys.argv[2]

## -------------------------------------------------------
## Scan input file and write output
## -------------------------------------------------------
fh_in  = open(vcf, 'r')
fh_out = open(out, 'w')

for line in fh_in.readlines():

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

	# Skip C -> T mutations for bisulfite transformation
	if ref == 'C' and alt_geno == 'T':
		continue

	# Skip G -> A mutations for bisulfite transformation in oposite strand
	if ref == 'G' and alt_geno == 'A':
		continue

	# Write line for position
	values = [chr, pos, ref, alt_geno]
	fh_out.write(' '.join(values))
	fh_out.write('\n')

## ----------------------------------------------------
## Close file handles
## ----------------------------------------------------
fh_in.close()
fh_out.close()
