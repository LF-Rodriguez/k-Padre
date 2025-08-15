## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 
##                              k'Padre-rna
##
##	idReads: This script parse a an alignment file containing reads
##		alignmed to a dipoloid genome and assignes genome or origin
##		based on a thresshold difference in alignment score.
##		
##		kPadre_RNA_idReads.py alignment.sam read_assignmen.txt
##
##		input: sam file sorted by bames
##		output read assignment table
##
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import sys
import copy
import subprocess
import re
import shutil


def timeReport():
	'''Returns date and time'''

	from datetime import datetime
	now = datetime.now()
	current_time = now.strftime("[%y/%m/%d %H:%M:%S]: ")
	return current_time

def is_read_mapped(flag):
	''' Returns True if the read is mapped, False if unmapped.'''
	return (flag & 0x4) == 0

def is_read_paired(flag):
	''' Returns True if the read is paired, Flase if unpaired. '''
	return (flag & 0x1) != 0

def is_read_firstPair(flag):
	''' Return True if the read is first in pair '''
	return (flag & 0x40) !=0

def is_read_mate_mapped(flag):
	''' Returns True if read mate is mapped ''' 
	return (flag & 0x8) == 0

def classifyAlignment(flag):
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
		
	if classifyAlignment(flag) == 'unpaired':

		for field in alignment_fields[11:]:

			if field.startswith("AS:i:"):
				AS = int(field.split(':')[2])
				break

		group = 'unpaired'
		YS    = None

	if classifyAlignment(flag) == 'unmated':
		for field in alignment_fields[11:]:

			if field.startswith("AS:i:"):
				AS = int(field.split(':')[2])
				break

		group = 'unmated'
		YS    = None

	if classifyAlignment(flag) == 'mated':
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
		alignment_info['pos'] = position

	except:
		alignment_info['group'] = 'problematic format'


	# Return results
	return alignment_info

def parseSam(filename):
	''' Parse sam file and collects relevant information 
	about each alignment. Returns dictionary with read 
	IDs as keys'''

	# create Sam file if needed 
	if filename.endswith('.bam'):
		NAME=re.sub(r'\.bam$', '', filename)
		CMD=f"samtools view {filename} > {NAME}.KPtmp"

		# Create SAM
		print(f"{timeReport()} Decompressing file with Samtools")
		subprocess.call(CMD, shell = True)
		filename=f'{NAME}.KPtmp'

	# Create file handle to write out problematic alignments
	NAME=re.sub(r'\.sam$', '', filename)
	out_handle = open(f'KP_{NAME}_problematic_alignments.txt', 'w')

	# Start dictionary and counters
	alignments = {}
	processed = 0
	failed = 0

	# Parse alignment file
	print(f"{timeReport()} Loading read alignment information")
	samfile = open(filename,'r')

	for line in samfile:
		if line.startswith('@'):
			continue

		fields = line.strip().split('\t')
		flag = int(fields[1])
		processed += 1

		# Skip reads that are not mapped ar that are second mate
		if not is_read_firstPair(flag):
			continue
		if not is_read_mapped(flag):
			continue

		# Collect alignment information
		alignment_info = collect_alignment_information(fields)

		# Discard problematic formated alignments
		if alignment_info['group'] == 'problematic format':

			failed += 1
			out_handle.write(f'Problem in read alignment format: {" ".join(fields)}')
			continue

		else:
			readID = alignment_info['id']
			genome = alignment_info['genome']
			AS     = alignment_info['AS']
			YS     = alignment_info['YS']
			group  = alignment_info['group']
			chrom  = alignment_info['chromosome']
			pos    = alignment_info['position']
			score  = AS + YS if group == 'mated' else AS
			summary = (genome, group, score, chrom, pos)

		# Store information if read exists
		if readID in alignments.keys():
			alignments[readID] = alignments[readID] + [summary]

		# Create new record for read in dictionary
		else:
			alignments[readID] = [summary]
	
	# close file handle 
	out_handle.close()

	return (alignments,processed,failed)

def AssignGenome(alignments, G1, G2):
	'''Compares alignment scores between genomes and assigns read
		To genome with higest score.

		alignments = dictionary with alignment information per read'''

	# Create copy of dictionary and start counters
	assignment = copy.deepcopy(alignments)
	nG1 = 0
	nG2 = 0
	nU  = 0

	# Loop through reads
	for read in alignments.keys():

		# Extract and group alignment information by genome
		readType = 'SE' if alignments[read][0][1] == 'unpaired' else 'PE'
		groups_G1 = [x[1] for x in alignments[read] if x[0] == G1]
		groups_G2 = [x[1] for x in alignments[read] if x[0] == G2]
		scores_G1 = [x[2] for x in alignments[read] if x[0] == G1]
		scores_G2 = [x[2] for x in alignments[read] if x[0] == G2]
		chrom_G1  = [x[3] for x in alignments[read] if x[0] == G1]
		chrom_G2  = [x[3] for x in alignments[read] if x[0] == G2]
		pos_G1    = [x[4] for x in alignments[read] if x[0] == G1]
		pos_G2    = [x[4] for x in alignments[read] if x[0] == G2]		

		# If read only aligned to one genome
		genomes = list(set([x[0] for x in alignments[read]]))
		
		if len(genomes) == 1:
			maxG1 = max(scores_G1) if len(scores_G1) > 0 else 'none'
			maxG2 = max(scores_G2) if len(scores_G2) > 0 else 'none'

			# Identify position of max socres 
			if maxG1 != 'none':
				idxG1 = maxG1.index(maxG1)
			elif maxG1 == 'none':
				idxG1 = 'none'

			if maxG2 != 'none':
				idxG2 = maxG2.index(maxG2)
			elif maxG1 == 'none':
				idxG2 = 'none'
			
			# Update counters and choose coordinate 
			if genomes[0] == G1:
				nG1  += 1
				idx  = maxG1.index(maxG1)
				pos  = pos_G1[idx]
				chro = chrom_G1[idx]

			if genomes[0] == G2:
				nG2  += 1
				idx  = max_G2.index(maxG2)
				pos  = pos_G2[idx]
				chro = chrom_G2[idx]

			# Create asignment
			assignment[read] = assignment[read] + [(genomes[0],maxG1, maxG2, chro, pos, 'UGM')]

		# If read mapped to two genomes
		if len(genomes) == 2:

			maxG1 = max(scores_G1)
			maxG2 = max(scores_G2)
		
			# ... and Read is Single end or only one end mapped in all cases:
			if readType == 'SE' or 'mated' not in groups_G1 + groups_G2:
			
				if maxG1 - maxG2 > 0:
					nG1+= 1
					idx  = maxG1.index(maxG1)
					pos  = pos_G1[idx]
					chro = chrom_G1[idx]
					assignment[read] = assignment[read] + [(G1,maxG1, maxG2, chro, pos, 'SEHS')]
					

				if maxG2 - maxG1 > 0:
					nG2+= 1
					idx  = maxG2.index(maxG2)
					pos  = pos_G2[idx]
					chro = chrom_G2[idx]
					assignment[read] = assignment[read] + [(G2,maxG1, maxG2, chro, pos,'SEHS')]
					

				else:
					assignment[read] = assignment[read] + [('undetermined',maxG1, maxG2,'na','na','ES')]
					nU += 1

			# ... or if mated alginment is present in only one genome:
			elif 'mated' in groups_G1 and 'mated' not in groups_G2:
				nG1+= 1
				idx   = maxG1.index(maxG1)
				pos   = pos_G1[idx]
				chro  = chrom_G1[idx]
				assignment[read] = assignment[read] + [(G1, maxG1, maxG2, chro, pos, 'UPA')]

			elif 'mated' in groups_G2 and 'mated' not in groups_G1:
				nG2+= 1
				idx  = maxG2.index(maxG2)
				pos  = pos_G2[idx]
				chro = chrom_G2[idx]
				assignment[read] = assignment[read] + [(G2, maxG1, maxG2, chro, pos, 'UPA')]


			# ... Lastly, if there are mated alignments in both
			elif 'mated' in groups_G1 and 'mated' in groups_G2:
				scores_G1_m = [scores_G1[i] for i in range(0,len(scores_G1)) if groups_G1[i] == 'mated']
				scores_G2_m = [scores_G2[i] for i in range(0,len(scores_G2)) if groups_G2[i] == 'mated']
			
				maxG1 = max(scores_G1_m)
				maxG2 = max(scores_G2_m)

				if maxG1 - maxG2 > 0:
					idx  = maxG1.index(maxG1)
					pos  = pos_G1[idx]
					chro = chrom_G1[idx]
					assignment[read] = assignment[read] + [(G1,maxG1, maxG2, chro, pos,'PEHS')]
					nG1+= 1 

				if maxG2 - maxG1 > 0:
					idx  = maxG2.index(maxG2)
					pos  = pos_G2[idx]
					chro = chrom_G2[idx]
					assignment[read] = assignment[read] + [(G2,maxG1, maxG2, chro, pos,'PEHS')]
					nG2+= 1

				else:
					assignment[read] = assignment[read] + [('undetermined',maxG1, maxG2,'NA','NA','PEES')]
					nU+= 1

	return (assignment, nG1, nG2, nU)

def writeTable(assignments, fname):
	'''Writes a text file with the results of the analysis of
	   the AssignGenome() function '''
 
	fh = open(f'KP_{fname}_readAssignment.txt', 'w')
 
	fh.write('ReadID\tmaxG1\tmaxG2\tgenome\tchr\tposition\tgroup\n')

	for k,v in assignments.items():
		read = k
		genome = v[-1][0]
		maxG1  = v[-1][1]
		maxG2  = v[-1][2]
		chrom  = v[-1][3]
		pos    = v[-1][4]
		group  = v[-1][5]
 
		fh.write(f'{read}\t{maxG1}\t{maxG2}\t{genome}\t{chrom}\t{pos}\t{group}\n')

	fh.close()

def main():

	# Catch wrong syntax
	if len (sys.argv) < 4:
		print ("Usage: kPadre_RNA_idReads.py input.sam prefix1 prefix2")
		sys.exit()
		
	# Check for dependencies 
	if shutil.which("samtools") is None:
    	sys.exit("Error: samtools is not installed or not in PATH. Please install it before running this script.")
	
	# obtain user arguments
	filename = sys.argv[1]
	G1 = sys.argv[2]
	G2 = sys.argv[3]
	fname=re.sub(r'\.[bs]am$', '', filename)

	print(f"{timeReport()} Processing file: {filename}")

	try:
		## Load read IDs and alignment metrics to dictionary
		alignments,p,f = parseSam(filename)

		print(f"{timeReport()} finished loading alignments")
		print(f"{timeReport()} -----Total alignments processed: {p}")
		print(f"{timeReport()} -----Reads discarded: {f}")
		print(f"{timeReport()} Starting genome assigenment")
		
		## Assign genome of origin to read
		assignments, nG1, nG2, nU = AssignGenome(alignments, G1, G2)
		
		print(f"{timeReport()} Finished genome assignment")
		print(f"{timeReport()} -----Total Reads assigned: {nG1 + nG2 + nU}")
		print(f"{timeReport()} --------Reads assigned to Genome 1: {nG1}")
		print(f"{timeReport()} --------Reads assigned to Genome 2: {nG2}")
		print(f"{timeReport()} --------Reads unassigned: {nU}")
		print(f"{timeReport()} Writing results to file")


		# Remove emporary files
		CMD='rm *.KPtmp'
		subprocess.call(CMD, shell = True)

		# Write results table
		writeTable(assignments, fname)

		print(f"{timeReport()} analsis completed")

	except FileNotFoundError:
		print(f"Error: File '{filename}' not found.")

if __name__ == "__main__":
	main()
