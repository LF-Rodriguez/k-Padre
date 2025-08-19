## -------------------------------------------------------
## QuePadre_v1.: Mark bisulfite transforemd reads in a sam
##		 file base on evidence of parental origin
##		 V1.1 breaks dictionary of informative SNPs
##		 in smaller volumes to speed up search
## -------------------------------------------------------

## -------------------------------------------------------
## Set up environment
## -------------------------------------------------------

# import modules .........................................
import sys
import re
import subprocess
import math

## -------------------------------------------------------
## Create funcitons
## -------------------------------------------------------

# Load table of snps - dividing dictionary by "Volumes" ...... 
def indexSNP_table(SNPs_file, vSize = 2000000):
	' Proces a file with reference and alternative SNP '
	' information to collect data as a list, and estimate '
	' de number of "volumes" to create an encyclopedia '
	' of SNPs '

	all_snps     = []
	volume_sizes = {}
	
	# Store raw information ................................................
	with open(SNPs_file, 'r') as fh:
		for line in fh.readlines():
			line = line.strip().split(' ')
			Chr  = line[0]
			pos  = int(line[1])
			G1   = line[2]
			G2   = line[3]
			vol  = math.floor(pos/vSize)
			
			# Store SNP in bank
			all_snps.append((Chr, pos, G1, G2, vol))
			
			# Update count of items by volume
			volume = Chr + '_' + str(vol)
			volume_sizes[volume] = volume_sizes.get(volume, 0) + 1
	
	# Return ................................................................
	nVol = len(volume_sizes)
	return (all_snps, volume_sizes, nVol, vSize)

def createChromopedia(all_snps, volume_sizes, nVol, margin = 50):
	'Uses the information from indexSNP_table to create'
	'an "enciclopedia" of SNPs, consisitng of multiple'
	'dictionaries with chr, pos and alleles divided by '
	'"Volumes" default to in volume every 1Mbp'
	''
	'chromopedia:['
	'	chr 1_0 <- {0010 : ["T","G"]}' 
	'	chr 1_1 <- {1025 : ["A","G"], 1045 : ["T","A"]}'
	'	chr 1_2 <- {2021 : ["A","G"], 2021 : ["T","A"]}'	
	
	chromopedia = {}
	cut = 0
	for vol in list(volume_sizes.keys()):
		name   = str(vol)
		size   = volume_sizes[vol]
		bottom = cut
		top    = (cut + size + margin) if cut == 0 else (cut + size + (2*margin))
		SNPset = all_snps[bottom:top]
		
		# Update cutting position
		cut = top - (2 * margin)
		
		# Create volume 
		chromopedia[name] = {}
		
		# Feed volume with SNPs
		for snp in SNPset:
			pos = snp[1]
			G1  = snp[2]
			G2  = snp[3]
			chromopedia[name][pos] = [G1, G2]
			
	return (chromopedia)

# define if read is mapped ...............................
def is_mapped(flag):
        'Returns True if read is mapped'

        binary = bin(flag)
        num    = binary[-3]
        if num == '0':
                return True
        else:
                return False

# Convert cigar string to ranges of coords ...............
def CigarToCoords(chr, start, cigar):
	'Returns tuple of bases covered by aligment'
	'(chr,[ R1, R2, R3 .. Rn],[G1,G2,G3, .. Gn])'

	CigarNumbers = re.split('[A-Z]', cigar)[:-1]
	CigarLetters = re.split('[0-9]+', cigar)[1:]

	GRunning = start
	RRunning = 0

	GCoords = []
	RCoords = []

	for i in range(0,len(CigarNumbers)):
		num  = int(CigarNumbers[i])
		type = CigarLetters[i]

		if type == 'M':
			GCoords += [x for x in range(GRunning, (GRunning + num))]
			RCoords += [x for x in range(RRunning, (RRunning + num))]
			GRunning += num
			RRunning += num

		if type == 'N':
			GRunning += num

		if type == 'D':
			GRunning += num

		if type == 'S':
			RRunning += num

		if type == 'I':
			RRunning += num
			
	return((chr, RCoords, GCoords))

# Retreive sites overlaping  diagnostic snp ..............
def OverDiagnostic(readSeq, coords, DiagSites, volumeSize = 2000000):
	'Return lists of tupples with diagnistic sites '
	'covered by alignment and read base calls covering '
	'diagnostic sites:'
	''
	'overDiagnostic = [(chr, pos, base),(chr,pos,base)]'

	Chr      = coords[0]
	RCoords  = coords[1]
	GCoords  = coords[2]
	sites    = list(DiagSites.keys())
	overLap  = list(set(GCoords) & set(sites))
	idx      = [GCoords.index(x) for x in overLap]

	overDiagnostic = []

	for i in idx:

		base = readSeq[RCoords[i]]
		pos  = GCoords[i]
		overDiagnostic.append((Chr, pos, base))
	
	return(overDiagnostic)

# Count snps suporting each genotype ......................
def SuportForAlleles(Overlap, DiagSites):
	'Returns a count of snps suporting each allele in '
	'alignment'

	Genome1 = 0
	Genome2 = 0

	for diagnostic in Overlap:
		Chr    = diagnostic[0]
		pos    = diagnostic[1]
		base   = diagnostic[2]
		
		alleles = DiagSites[pos]
		G1      = alleles[0]
		G2      = alleles[1]

		if base == G1:
			Genome1+=1
		if base == G2:
			Genome2+=1

	return [Genome1, Genome2]

# Generate time reports ..................................
def timeReport():
	'Returns date and time'

	from datetime import datetime
	now = datetime.now()
	current_time = now.strftime("[%y/%m/%d %H:%M:%S]: ")
	return current_time

## -------------------------------------------------------
## Read input files and collect file information
## -------------------------------------------------------

# read user arguments ....................................
inBam = sys.argv[1]
inKey = sys.argv[2]
format = sys.argv[3]

BamfileName = inBam.split('/')[-1]
keyfileName = inKey.split('/')[-1]

# Grab file id  ..........................................
fileID = inBam[:-4]

## -------------------------------------------------------
## Generate indexed enciclopedia of SNPs
## -------------------------------------------------------

# Update progres .........................................
print( ' -----------------')
print("K'Padre run")
print('Alignment file: ' + BamfileName )
print('Variants data base: ' + keyfileName )
print(' ')
print(timeReport() + 'Generating indexed data base of varinats')

# Index SNP database .....................................
SNPindex   = indexSNP_table(inKey)
allSNPS    = SNPindex[0]
BlockSizes = SNPindex[1]
nVolumes   = SNPindex[2]
volumeSize = SNPindex[3]

# Create enciplopedia of SNPs ............................
chromoPedia = createChromopedia(allSNPS, BlockSizes, nVolumes)

# Update progres .........................................
print(timeReport() + 'Finished generating data base')

## -------------------------------------------------------
## Write temorary files
## -------------------------------------------------------

# Update progres .........................................
print(timeReport() + 'Writing temporary files')

# Print header ...........................................
comand_01 = 'samtools view -H ' + inBam + ' > ' + fileID + '_KP.tmp01.header'
subprocess.call(comand_01, shell = True)

# Print original alignments ..............................
comand_02 = 'samtools view ' + inBam + ' > ' + fileID + '_KP.tmp02.sam'
subprocess.call(comand_02, shell = True)

# Create file handles ....................................
fh1 = open(fileID + '_KP.tmp02.sam', 'r')
fh2 = open(fileID + '_KP.tmp03_NewAlignments.sam', 'w')

# Update progres .........................................
print(timeReport() + 'Procesing alignments')

# Start counters ........................................
records = 0
assigned = 0
parent1 = 0
parent2 = 0

## -------------------------------------------------------
## Run K' Padre algorithm
## -------------------------------------------------------

# Process alignments .....................................
for line in fh1.readlines():

	# Update counter
	records += 1

 	# Grab alignment information
	line     = line.strip().split()
	read     = line[0]
	flag     = int(line[1])
	chr      = line[2]
	start    = int(line[3])
	cigar    = line[5]
	sequence = line[9]
	block    = math.floor( start / volumeSize)
	keyName  = str(chr) + '_' + str(block)

	# Process alignment
	if keyName in chromoPedia.keys():
		if is_mapped(flag):
			snps      = chromoPedia[keyName]
			coords    = CigarToCoords(chr, start, cigar)
			diagBases = OverDiagnostic(sequence, coords, snps, volumeSize)
			support   = SuportForAlleles(diagBases, snps)

			# Choose parent
			if support[0] > support[1]:
				parent    = 1
				parent1  += 1
				assigned += 1

			if support[0] < support[1]:
				parent    = 2
				parent2  += 1
				assigned += 1

			if support[0] == support[1]:
				parent = 0

			# Collect updated info 
			po = 'po:i:' + str(parent)
			ns = 'ns:i:' + str(support[0] + support[1])
			g1 = 'g1:i:' + str(support[0])
			g2 = 'g2:i:' + str(support[1])
			old = '\t'.join(line)

		# Alignment does not overlap any block of SNPs
		else:
			po  = 'po:i:0'
			ns  = 'ns:i:0'
			g1  = 'g1:i:0'
			g2  = 'g2:1:0'
			old = '\t'.join(line)
		
	
		# Write output
		fh2.write(old + '\t' + ns + '\t' + g1 + '\t' + g2 + '\t' + po + '\n')

# Update progres .........................................
print(timeReport() + 'Finished processing alignments')
print(timeReport() + 'Generating final alignment file')
			
# Close file handles
fh1.close()
fh2.close()

# Concatenate header end new alignments ...................
comand_03 = 'cat ' + fileID + '_KP.tmp01.header ' + fileID + '_KP.tmp03_NewAlignments.sam > ' + fileID + '_KP.sam'
subprocess.call(comand_03, shell = True)

# Create Bam if requested ................................
if format == 'bam':
	comand_04 = 'samtools view -b ' + fileID + '_KP.sam > ' + fileID + '_KP.bam'
	comand_05 = 'rm ' + fileID + '_KP.tmp* ' + fileID + '_KP.sam'

	subprocess.call(comand_04, shell = True)
	subprocess.call(comand_05, shell = True)

# Remove temporary files ..................................
if format == 'sam':
	comand_04 = 'rm ' + fileID + '_KP.tmp* '
	subprocess.call(comand_04, shell = True)

# Write final summary ....................................
print(timeReport() + 'Finished run \n')
print('Summary report: \n')
print('--------------- \n \n')
print('Alignments processed: '+ str(records))
print('Alignments assigned to parents: ' + str(assigned))
print('\t Assigned to parent 1: ' + str(parent1))
print('\t Assigned to parent 2: ' + str(parent2))

