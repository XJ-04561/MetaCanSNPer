'''
ParseXMFA is a script that parse xmfa files and extract all SNPs
	The script allows an option using a flanking score that limits
	SNPs near edges of a block to have a impact on classifying decisions.
'''

#__name__="ParseXMFA"
import os
import re
import mmap
import logging
logger = logging.getLogger(__name__)

reference_genomes = ["FSC200","SCHUS4.1","SCHUS4.2","OSU18","LVS","FTNF002-00"]

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
		return -1

	def find(self, pos : int):
		'''Return where in the associated file object that this position corresponds to.'''
		return self[pos]
	


class ParseXMFA(object):
	entryHeader : re.Pattern = re.compile(b"^[>] [0-9]+:(?P<start>[0-9]+)-(?P<stop>[0-9]+) - .*", re.MULTILINE)
	filename : str
	_mmap : mmap.mmap
	coverage : Map

	def __init__(self, filename : str=None, mode : str="r"):
		if filename is not None:
			self.open(filename, mode)
			self.mapBlocks()
	
	def __contains__(self, item):
		return item in self.coverage
	
	def __getitem__(self, key : list[int]|int):
		offset = self.index(key)
		if type(offset) is list:
			return [self._mmap[o] for o in offset]
		else:
			return self._mmap[offset]
				
	def open(self, filename, mode : str):
		self.filename = filename
		self.file = open(self.filename, "r+b") # Makes the file writeable, but we do not intend to edit it
		self.file.readline()
		self.coverage = Map(self.file.newlines)
		self.file.seek(0)

		if mode == "rw":
			self._mmap = mmap.mmap(self.file.fileno(), 0)
		elif mode == "r":
			self._mmap = mmap.mmap(self.file.fileno(), 0, prot=mmap.PROT_WRITE)
		elif mode == "w":
			self._mmap = mmap.mmap(self.file.fileno(), 0, prot=mmap.PROT_READ)
		else:
			raise ValueError("'{}' is not a valid mode for opening files.".format(mode))
	
	def mapBlocks(self):
		for m in self.entryHeader.finditer():
			filepos = m.end()+2
			ntRange = range(int(m.group("start")), int(m.group("stop")))
			self.coverage += filepos, ntRange
		
	def index(self, pos : int, *mPos : int):
		'''Takes more than 0 integer positions in the sequence covered by the XMFA and returns the byte position(s) in
		the XMFA file that correspond to the given position in the aligned sequence. If only one position is given the
		return type is int. If more than one positions are ggiven, a list of byte positions is returned.'''
		if len(mPos) == 0:
			return self.coverage[pos]
		else:
			return [self.coverage[p] for p in [pos]+mPos]

if __name__=="__main__":
	parser = ParseXMFA(os.path.join("CanSNPer2", "test.xmfa"))
	