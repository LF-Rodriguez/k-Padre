#!/usr/bin/env python3

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 
##                              kPadre twoGenomes 
##
##		This script parse a an alignment file containing reads
##		alignmed to a dipoloid genome and assigns genome of origin
##		based on a thresshold difference in alignment score.		
##
##		input: alignment to diploid genome (SAM / BAM)
##		output: read assignment table, problematec reads table
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Import 
import os
import sys
import copy
import subprocess
import re
import shutil
import logging
import os
import glob

# Add project root to Python path to import custom modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.sam_utils import *
from utils.io_utils import *

def classify_alignment(flag):
	''' Clasifies alignment in three categories based on mate presence
		and alignment
	
		unpaired = does not have a mate
		mated    = mate exist and was mapped
		unmated  = mate exists and was not mapped
	'''

	if not is_read_paired(flag):
		return 'unpaired'

	if is_read_mate_mapped(flag):
		return 'mated'

	else:
		return 'unmated'

def collect_alignment_information(alignment_fields):
	''' Reads in an alignment (line from SAM file), packed as a list of
		elements and collect all relevant information about alignment.
		returns None if read is not aligned or dictionary of stats
		if it is aligned

		alignment_fields = list 
	'''

	alignment_info = {}

	read_id = alignment_fields[0]
	flag    = int(alignment_fields[1])
	genome  = alignment_fields[2][:3]
	chrom   = alignment_fields[2]
	position = alignment_fields[3] 
		
	if classify_alignment(flag) == 'unpaired':

		for field in alignment_fields[11:]:

			if field.startswith("AS:i:"):
				AS = int(field.split(':')[2])
				break

		group = 'unpaired'
		YS    = None

	if classify_alignment(flag) == 'unmated':
		for field in alignment_fields[11:]:

			if field.startswith("AS:i:"):
				AS = int(field.split(':')[2])
				break

		group = 'unmated'
		YS    = None

	if classify_alignment(flag) == 'mated':
		for field in alignment_fields[11:]:

			if field.startswith("AS:i:"):
				AS = int(field.split(':')[2])
				
			if field.startswith("YS:i:"):
				YS = int(field.split(':')[2])

		group = 'mated'

	# Store alignment information
	try:

		alignment_info['id'] = read_id
		alignment_info['genome'] = genome
		alignment_info['AS'] = AS
		alignment_info['YS'] = YS
		alignment_info['group'] = group
		alignment_info['chromosome'] = chrom
		alignment_info['position'] = position

	except:
		alignment_info['group'] = 'problematic format'


	# Return results
	return alignment_info

def parse_sam(filename):
	'''Parse SAM or BAM file and collect relevant alignment information.
    	Returns a dictionary with read IDs as keys, along with counters.
    
    	Parameters:
    	- filename: Path to input .sam or .bam file
    	
    	Returns:
    	- dict: alignments keyed by read ID
    	- int: number of processed reads
    	- int: number of failed/problematic reads
    '''
    # Import libraries 
	import os
	import re
	import subprocess
	import logging

	# Start counters 
	logging.basicConfig(level=logging.INFO, format='%(message)s')
	alignments = {}
	processed = 0
	failed = 0
    
	# Fetch names and create output file names 
	base = os.path.basename(filename)
	name = re.sub(r'\.bam$', '', base)
	alg_file = f"{name}_alignments.KPtmp"
	prm_file = f"{name}_problematic_alignments.txt"
	
	# Create file of raw alignments, skiping unmapped reads
	logging.info(f"{time_report()} Decompressing file with samtools")

	with open(alg_file ,'w') as fh_out:
		subprocess.run(["samtools", "view", "-F", "0x4", filename], stdout=fh_out)
	
	# Parse Sam file and report problematic alignments
	with open(alg_file, 'r') as fh_in, open(prm_file, 'w') as fh_out:

		logging.info(f"{time_report()} Loading read alignment information")
		
		for line in fh_in:
			if line.startswith('@'):
				continue

			fields = line.strip().split('\t')

			# Exclude incomplete alignments
			if len(fields) < 11: # minimum number of SAM fields
				failed += 1
				processed += 1
				fh_out.write('Incomplete alignment: {" ".join(fields)}\n')
				continue

			# Obtain alignment flag
			flag = int(fields[1])
			processed += 1

			if is_read_paired(flag) and not is_read_firstPair(flag):
				continue

			# Collect alignment information
			alignment_info = collect_alignment_information(fields)

			# Discard problematic formated alignments
			if alignment_info['group'] == 'problematic format':
				failed += 1
				fh_out.write(f'Problem in read alignment format: {" ".join(fields)}\n')
				continue

			# Extract alignment information
			readID = alignment_info['id']
			genome = alignment_info['genome']
			AS = alignment_info['AS']
			YS = alignment_info['YS']
			group = alignment_info['group']
			chrom = alignment_info['chromosome']
			pos = alignment_info['position']
			score = AS + YS if group == 'mated' else AS
			summary = (genome, group, score, chrom, pos)

			# Add entry to dictionary of alignmens
			alignments.setdefault(readID, []).append(summary)

	return (alignments, processed, failed)

def assign_genome(alignments, G1, G2):
    """
    Compares alignment scores between genomes and assigns each read
    to the genome with the highest score.

    Parameters:
        alignments (dict): Dictionary with alignment information per read.
        G1, G2 (str): Genome identifiers.

    Returns:
        tuple: (assignment dict, count G1, count G2, count undetermined)
    """

    # Create copy of dictionary and start counters
    assignment = copy.deepcopy(alignments)
    nG1 = 0
    nG2 = 0
    nU = 0

    # Loop through reads
    for read in alignments.keys():

        # Extract alignment information by genome
        read_type = 'SE' if alignments[read][0][1] == 'unpaired' else 'PE'
        groups_G1 = [x[1] for x in alignments[read] if x[0] == G1]
        groups_G2 = [x[1] for x in alignments[read] if x[0] == G2]
        scores_G1 = [x[2] for x in alignments[read] if x[0] == G1]
        scores_G2 = [x[2] for x in alignments[read] if x[0] == G2]
        chrom_G1 = [x[3] for x in alignments[read] if x[0] == G1]
        chrom_G2 = [x[3] for x in alignments[read] if x[0] == G2]
        positions_G1 = [x[4] for x in alignments[read] if x[0] == G1]
        positions_G2 = [x[4] for x in alignments[read] if x[0] == G2]

        # Case 1: read only aligned to one genome
        genomes = list(set([x[0] for x in alignments[read]]))

        if len(genomes) == 1:
            maxG1 = max(scores_G1) if scores_G1 else 'none'
            maxG2 = max(scores_G2) if scores_G2 else 'none'

            # Identify position of max scores
            idxG1 = scores_G1.index(maxG1) if maxG1 != 'none' else 0
            idxG2 = scores_G2.index(maxG2) if maxG2 != 'none' else 0

            # Update counters and choose coordinate
            if genomes[0] == G1:
                nG1 += 1
                idx = scores_G1.index(maxG1)
                posG1 = positions_G1[idx]
                chrG1 = chrom_G1[idx]
                posG2 = 0
                chrG2 = 'none'

            if genomes[0] == G2:
                nG2 += 1
                idx = scores_G2.index(maxG2)
                posG2 = positions_G2[idx]
                chrG2 = chrom_G2[idx]
                posG1 = 0
                chrG1 = 'none'

            # Create assignment
            assignment[read] += [(genomes[0], maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'UGM')]

        # Case 2: Read mapped to two genomes
        if len(genomes) == 2:

            # Obtain maximum mapping scores
            maxG1 = max(scores_G1)
            maxG2 = max(scores_G2)

            # Collect mapping info from genome 1
            idxG1 = scores_G1.index(maxG1)
            posG1 = positions_G1[idxG1]
            chrG1 = chrom_G1[idxG1]

            # Collect mapping info for genome 2
            idxG2 = scores_G2.index(maxG2)
            posG2 = positions_G2[idxG2]
            chrG2 = chrom_G2[idxG2]

            # Sub-case 2-1: Read is single end or only one end mapped in both genomes
            if read_type == 'SE' or ('mated' not in groups_G1 + groups_G2):

                if maxG1 > maxG2:
                    nG1 += 1
                    assignment[read] += [(G1, maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'SEHS')]

                elif maxG2 > maxG1:
                    nG2 += 1
                    assignment[read] += [(G2, maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'SEHS')]

                else:
                    assignment[read] += [('undetermined', maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'SEES')]
                    nU += 1

            # Sub-case 2-2: Mated alignment is present in only one genome
            elif 'mated' in groups_G1 and 'mated' not in groups_G2:
                nG1 += 1
                assignment[read] += [(G1, maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'UPA')]

            elif 'mated' in groups_G2 and 'mated' not in groups_G1:
                nG2 += 1
                assignment[read] += [(G2, maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'UPA')]

            # Sub-case 2-3: There are mated alignments in both genomes
            elif 'mated' in groups_G1 and 'mated' in groups_G2:

                # Subset alignment info to include only scores of mated alignments
                scores_G1_m = [scores_G1[i] for i in range(len(scores_G1)) if groups_G1[i] == 'mated']
                scores_G2_m = [scores_G2[i] for i in range(len(scores_G2)) if groups_G2[i] == 'mated']
                chrom_G1_m = [chrom_G1[i] for i in range(len(chrom_G1)) if groups_G1[i] == 'mated']
                chrom_G2_m = [chrom_G2[i] for i in range(len(chrom_G2)) if groups_G2[i] == 'mated']
                positions_G1_m = [positions_G1[i] for i in range(len(positions_G1)) if groups_G1[i] == 'mated']
                positions_G2_m = [positions_G2[i] for i in range(len(positions_G2)) if groups_G2[i] == 'mated']

                # Select maximum scores
                maxG1 = max(scores_G1_m)
                maxG2 = max(scores_G2_m)

                # Collect mapping info from genome 1
                idxG1 = scores_G1_m.index(maxG1)
                posG1 = positions_G1_m[idxG1]
                chrG1 = chrom_G1_m[idxG1]

                # Collect mapping info for genome 2
                idxG2 = scores_G2_m.index(maxG2)
                posG2 = positions_G2_m[idxG2]
                chrG2 = chrom_G2_m[idxG2]

                if maxG1 > maxG2:
                    assignment[read] += [(G1, maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'PEHS')]
                    nG1 += 1

                elif maxG2 > maxG1:
                    assignment[read] += [(G2, maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'PEHS')]
                    nG2 += 1

                else:
                    assignment[read] += [('undetermined', maxG1, maxG2, chrG1, posG1, chrG2, posG2, 'PEES')]
                    nU += 1

    return assignment, nG1, nG2, nU

def main( filename, G1, G2):

	# Set loggin parameters
	logging.basicConfig(level=logging.INFO, format='%(message)s')

	# Check for dependencies 
	if shutil.which("samtools") is None:
		sys.exit("Error: samtools is not installed or not in PATH. Please install it before running this script.")

	base = os.path.basename(filename)
	fname=re.sub(r'\.[bs]am$', '', base)

	logging.info(f"{time_report()} running kPadre twoGenomes {filename} {G1} {G2}")
	logging.info(f"{time_report()} Processing file: {filename}")

	try:
		## Load read IDs and alignment metrics to dictionary
		alignments,p,f = parse_sam(filename)

		logging.info(f"{time_report()} finished loading alignments")
		logging.info(f"{time_report()} -----Total alignments processed: {p}")
		logging.info(f"{time_report()} -----Reads discarded: {f}")
		logging.info(f"{time_report()} Starting genome assigenment")
		
		## Assign genome of origin to read
		assignments, nG1, nG2, nU = assign_genome(alignments, G1, G2)
		
		logging.info(f"{time_report()} Finished genome assignment")
		logging.info(f"{time_report()} -----Total Reads assigned: {nG1 + nG2 + nU}")
		logging.info(f"{time_report()} --------Reads assigned to Genome 1: {nG1}")
		logging.info(f"{time_report()} --------Reads assigned to Genome 2: {nG2}")
		logging.info(f"{time_report()} --------Reads unassigned: {nU}")
		logging.info(f"{time_report()} Writing results to file")

		# Remove emporary files
		for temp_file in glob.glob("./*.KPtmp"):
			os.remove(temp_file)
		
		# Write results table
		write_scores_table(assignments, fname)

		logging.info(f"{time_report()} analsis completed")

	except FileNotFoundError:
		logging.info(f"Error: File '{filename}' not found.")

# Function for dispatcher kPadre.py
def run(filename, G1, G2):
	main(filename, G1, G2)

# Run as stand alone
if __name__ == "__main__":

	# Catch wrong syntax
	if len (sys.argv) < 4:
		print ("Usage: kPadre_twoGenomes input.sam prefix1 prefix2")
		sys.exit()

	# Obtain positional arguments
	filename = sys.argv[1]
	G1 = sys.argv[2]
	G2 = sys.argv[3]
   	
	main(filename, G1, G2)
