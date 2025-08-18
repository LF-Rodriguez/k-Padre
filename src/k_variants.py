## ------------------------------------------------------------------------------
## K-variants: Scan vcf of filtered genotypes to create
##	   a table of diagnostic SNPs for allele assignment of
##	   Whole Genome Bisulphite Seuencing (WGBS) data.
##
## ------------------------------------------------------------------------------

# Set environment ...............................................................

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__) , '..')))
import sys
import re
import subprocess
from kpadre.analysis import select_WGBS_Variants

# Define main function ..........................................................

def main():
	
	# Load vcf 
	vcf = inputVCF

	# Select variants
	variants = select_WGBS_Variants(vcf)

	# Save reults

# Entry point ....................................................................

if __name__ == '__main__':
	main()



