## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 
## 	Utilities for classifying sam alignments
## 
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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