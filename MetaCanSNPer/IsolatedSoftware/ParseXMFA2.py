'''
ParseXMFA is a script that parse xmfa files and extract all SNPs
	The script allows an option using a flanking score that limits
	SNPs near edges of a block to have a impact on classifying decisions.
'''

#__name__="ParseXMFA"
import re, mmap, logging

import MetaCanSNPer.modules.LogKeeper as LogKeeper
import MetaCanSNPer.modules.VCFhandler as openVCF

LOGGER = LogKeeper.getLogger(__name__)
LINEWIDTH = 80

class PseudoRange:
	start : int
	stop : int
	def __init__(self, start, stop):
		self.start = start
		self.stop = stop
	
	def __lt__(self, n):
		return self.start < n
	
	def __gt__(self, n):
		return self.stop > n

class Map:
	_ranges : list[range]
	_filepos : list[int]
	_newline : str

	def __init__(self, newline):
		self._ranges = []
		self._filepos = []
		self._newline = newline

	def __contains__(self, pos):
		for r in self._ranges:
			if pos >= r.start:
				if pos <= r.stop:
					return True
				break
		return False

	def __iadd__(self, newRange : tuple[int,range]):
		for i, r in enumerate(self._ranges):
			# Overlaps should not happen?
			# if newRange.stop > r.start:
			# 	if newRange.start < r.stop:
			# 		# Overlapping
			# 		self._ranges[i] = range(min(newRange.start, r.start), max(newRange.stop, r.stop))
			if newRange[1].start > r.start:
				self._filepos.insert(i, newRange[0])
				self._ranges.insert(i, newRange[1])
	
	def __getitem__(self, pos : int):
		'''Return where in the associated file object that this position corresponds to.'''
		for i, r in enumerate(self._ranges):
			if pos > r.start:
				if pos > r.stop:
					break
				else:
					return self._filepos[i] + pos - r.start + len(self._newline)*(pos - r.start)//LINEWIDTH
		return "-"

	def find(self, pos : int):
		'''Return where in the associated file object that this position corresponds to.'''
		return self[pos]
	


class ParseXMFA:
	entryHeader : re.Pattern = re.compile(b"^[>] (?P<seqID>[0-9]+):(?P<start>[0-9]+)-(?P<stop>[0-9]+) - .*", re.MULTILINE)
	reference : int
	filename : str
	_mmap : mmap.mmap
	coverage : list[Map]

	def __init__(self, filename : str=None, mode : str="r", referenceID : int=1):
		if filename is not None:
			self.open(filename, mode, referenceID)
			self.mapBlocks()
	
	def __contains__(self, item):
		return any(item in map for map in self.coverage)
	
	def __getitem__(self, key : list[int]|int):
		offset = self.index(key)
		if type(offset) is list:
			return map(self._mmap.__getitem__, offset)
		else:
			return self._mmap[offset]
				
	def open(self, filename, mode : str, referenceID : int):
		self.filename = filename
		self.referenceID = referenceID
		LOGGER.debug("open('{}', 'r+b')".format(filename))
		self.file = open(self.filename, "r+b") # Makes the file writeable, but we do not intend to edit it
		self.file.readline()
		self.newlines = self.file.newlines
		self.coverage = {}
		self.file.seek(0)

		if mode == "rw":
			self._mmap = mmap.mmap(self.file.fileno(), 0)
		elif mode == "r":
			self._mmap = mmap.mmap(self.file.fileno(), 0, prot=mmap.PROT_WRITE)
		elif mode == "w":
			self._mmap = mmap.mmap(self.file.fileno(), 0, prot=mmap.PROT_READ)
		else:
			raise ValueError(f"{mode!r} is not a valid mode for opening files.")
	
	def mapBlocks(self):
		'''Can handle multiple alignments'''
		block = {}
		for m in self.entryHeader.finditer():
			ntRange = range(int(m.group("start")), int(m.group("stop")))
			filepos = m.end()+len(self.newlines)
			if int(m.group("seqID")) in block:
				if self.referenceID not in block:
					LOGGER.error("XMFA file must have blocks that match the given referenceID sequence.")
					raise ValueError("XMFA file must have blocks that match the given referenceID sequence.")
				for i in block:
					if i != self.referenceID:
						if i not in self.coverage: self.coverage[i] = Map(self.newlines)
						self.coverage[i] += block[i][0], block[self.referenceID][1]
						
				self.block = {}
			else:
				block[int(m.group("seqID"))] = (filepos, ntRange)
	
	def index(self, pos : int, *mPos : int):
		'''Takes more than 0 integer positions in the sequence covered by the XMFA and returns the byte position(s) in
		the XMFA file that correspond to the given position in the aligned sequence. If only one position is given the
		return type is int. If more than one positions are ggiven, a list of byte positions is returned.'''
		if len(mPos) == 0:
			return [map[pos] for map in self.coverage]
		else:
			return [[map[p] for map in self.coverage] for p in [pos]+mPos]

if __name__=="__main__":
	import argparse
	aparse = argparse.ArgumentParser()
	aparse.add_argument("XMFAfile", metavar="XMFAfile", help="Alignment .XMFA file to load.")
	aparse.add_argument("SNPfile", metavar="SNPfile", help="A .vcf file of the SNPs to be called.")
	aparse.add_argument("-o", "--output", metavar="outFile", help="Name of the .vcf file to be created.")
	aparse.add_argument("-rID", "--referenceID", metavar="referenceID", help="The reference sequence used to create the .XMFA file.")

	args = aparse.parse_args()

	print("Creating Index.")
	parser = ParseXMFA(filename=args.XMFAfile, referenceID=args.referenceID)

	print(f"Opening {args.SNPfile!r}.")
	snpFile = openVCF(args.SNPfile, "r")
	
	print(f"Creating {args.outFile!r}.")
	outFile = openVCF(args, "w", referenceFile="")
	for entry in snpFile:
		outFile.add(POS=entry.POS, REF=parser[entry.POS])
	outFile.close()
	print("Done!")