#!/usr/bin/env python3

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 
##                              			kPadre
##
##		kPadre provides utilities for assigning allelic origin of sequencing reads. Allelic asignment 
##		can be performed thrugh two approaches:
##
##		- Single reference: using one reference genome assembly combined with a set of diagnostic 
##		  snps to distinguish parental lines
##
##		- Double reference: using a diploid genome assembly containing sequences from both parents.
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# import protocols 
import sys
import argparse
import os
from protocols import kPadre_twoGenomes, kPadre_selectSNPs, kPadre_createGenome

# Add current location of master script
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Define parser function
def getParsers():
	parser = argparse.ArgumentParser(description = "Asign allelic origin to reads using k-Padre")

	subparsers = parser.add_subparsers(
		title = "Protocols", dest = "mode", required = True, help = 'Choose protocol to use',
	)

	# Parser for single-reference protocol
	one_parser = subparsers.add_parser("oneGenome", help = "Use one genome reference and diagnostic SNPs")
	one_parser.add_argument("-a", required = True, metavar = "alignment.bam", help = "Path to genome alignment file")
	one_parser.add_argument("-s", required = True, metavar = "snps.txt", help = "Path to diagnostic SNPs file" )

	# Parser for two-genomes protocol
	two_parser = subparsers.add_parser("twoGenomes", help = "Use alignment to a diploid genome")
	two_parser.add_argument("-a", required = True, metavar = "alignment.bam", help = "path to diploid genome alignment file")
	two_parser.add_argument("-p1", required = True, metavar = "G1", help = "prefix for parent 1")
	two_parser.add_argument("-p2", required = True, metavar = "G2", help = "prefix for parent 2")

	# Parser for generating diploid genome
	genome_parser = subparsers.add_parser("createGenome", help = "Generate diploid genome file")
	genome_parser.add_argument("-g1", required = True, metavar = "genome1.fasta", help = "Path to genome of parent 1")
	genome_parser.add_argument("-g2", required = True, metavar = "genome2.fasta", help = "Path to genome of parent 2")
	genome_parser.add_argument("-p1", required = True, metavar = "G1", help = "prefix for parent 1")
	genome_parser.add_argument("-p2", required = True, metavar = "G2", help = "prefix for parent 2")

	# Parser for generating SNP file
	genome_parser = subparsers.add_parser("selectSNPs", help = "Generate SNPs file")
	genome_parser.add_argument("-v", required = True, metavar = "variants.vcf", help = "Path to VCF file with variants")
	genome_parser.add_argument("-bs", required = False, action="store_true",
		help = "include if processing bisulphite sequence data"
	)
	
	return parser

# Define main function
def main():
	parser = getParsers()

	# Check if help is  needed
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(0)

	args = parser.parse_args()

	if args.mode == "oneGenomes":
		kPadre_oneGenome.run(args.a, args.s)

	elif args.mode == "twoGenomes":
		kPadre_twoGenomes.run(args.a, args.p1, args.p2)

	elif args.mode == "createGenome":
		kPadre_createGenome.run(args.g1, args.g2, args.p1, args.p2)

	elif args.mode == "selectSNPs":
		kPadre_selectSNPs.run(args.v, args.bs)

	else:
		sys.exit(0)

if __name__ == "__main__":
	main()
