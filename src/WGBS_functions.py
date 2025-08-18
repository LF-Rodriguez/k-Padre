## ------------------------------------------------------------------------------
## K-variants: Funcitons for allelic assignment using diagnostic SNPs
##
## ------------------------------------------------------------------------------

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

def is_mapped(flag):
        'Returns True if read is mapped'

        binary = bin(flag)
        num    = binary[-3]
        if num == '0':
                return True
        else:
                return False

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

def timeReport():
	'Returns date and time'

	from datetime import datetime
	now = datetime.now()
	current_time = now.strftime("[%y/%m/%d %H:%M:%S]: ")
	return current_time