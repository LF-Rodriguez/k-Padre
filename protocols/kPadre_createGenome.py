#!/usr/bin/env python3

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 
##                              kPadre createGenome 
##
##		This script creates a diploid genome for kPadre.py twoGenomes
##		adding a prefix to each scaffold Id marking the origin before
##		merghing
##
##		input: genome1.fasta, genome2.fasta
##		output: siploid_genome.fasta
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import os
import sys
import subprocess

def createGenome(file1, file2, p1, p2):
	'''Concatenate two genomes adding a prefix to track
		origin
    
    	Parameters:
    	- file1: Path to genome 1
    	- file 2: File to genome 2
    	- p1 prefix for genome 1
    	- p2 prefix for genome 2 
    	
    	Prints:
    	- diploid genome to standar output
    '''

    # Add prefix to genome 1
	with open('g1.KPtmp', 'w') as fh_out:
		subprocess.run(["sed",f"s/>/>{p1}/g", file1], stdout = fh_out)

    # Add prefix to genome 2
	with open('g2.KPtmp', 'w') as fh_out:
		subprocess.run(["sed",f"s/>/>{p2}/g", file2], stdout = fh_out)

    # concatenate genomes
	subprocess.run(["cat", "g1.KPtmp",'g2.KPtmp'])

    # Delete temporary files
	subprocess.run(["rm", "g1.KPtmp","g2.KPtmp"])

# Define main function
def main(file1, file2, p1, p2):
	createGenome(file1, file2, p1, p2)

# Function for dispatcher kPadre.py
def run(file1, file2, p1, p2):
	main(file1, file2, p1, p2)

# Run as stand alone
if __name__ == "__main__":

	# Catch wrong syntax
	if len (sys.argv) < 5:
		print ("Usage: kPadre_createGenome genome1.fasta genome2.fasta prefix1 prefix2 > diploid_genome.fasta")
		sys.exit()

	# Obtain positional arguments
	file1 = sys.argv[1]
	file2 = sys.argv[2]
	p1 = sys.argv[3]
	p2 = sys.argv[4]
   	
	main(file1, file2, p1, p2)
