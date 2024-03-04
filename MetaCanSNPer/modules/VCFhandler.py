
from typing import TextIO, overload
import time, logging
try:
	import MetaCanSNPer.modules.LogKeeper as LogKeeper
except:
	import LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)


EXAMPLE_VCF_HEADER = """##fileformat=VCFv4.3
##fileDate={dateYYYYMMDD}
##source=MetaCanSNPer
##reference={referenceFile}
##contig=<ID={ID},length={length},assembly={assembly},md5={md5},species="{species}",taxonomy="{taxonomy}">
##phasing={phasing}
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=BQ,Number=1,Type=Integer,Description="RMS base quality at this position">
##FORMAT=<ID=MQ,Number=2,Type=Integer,Description="RMS mapping quality, e.g. MQ=52">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
"""

VCF_HEADER = """##fileformat=VCFv4.3
##fileDate={dateYYYYMMDD}
##source=MetaCanSNPer
##reference=file://{refPath}
##contig=<ID={chrom},URL=file://{refPath}>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""

COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
OPTIONAL_COLUMNS = ["FORMAT"]
VCF_ROW =       "{CHROM}	{POS}	{ID}	{REF}	{ALT}	{QUAL}	{FILTER}	{INFO}"
SEPARATORS = [(0, ";"), (1, ":"), (2, "|"), (3, ",")]

class RowDict: pass

def interpretIterable(string : str, seps=SEPARATORS):
	if string.startswith("[") and string.endswith("]"):
		string = string[1:-1]
		for i, s in seps:
			if s in string:
				return [interpret(subString, seps=seps[i+1:]) for subString in string.split(s)]
		return [interpret(string, seps=[])]
	else:
		for i, s in seps:
			if s in string:
				return [interpret(subString, seps=seps[i+1:]) for subString in string.split(s)]
		return string

def interpret(string : str|bytes, seps=SEPARATORS):
	"""Interprets anything with decimals in it or surrounded by [ ] as a list of items. Attempts to convert to int,
	then float, and then to decode with utf-8."""
	if type(string) is bytes:
		string : str = string.decode("utf-8")
	out = interpretIterable(string, seps=seps)
	if out != string:
		return out
	
	try:
		return int(string)
	except:
		pass

	try:
		return float(string)
	except:
		pass

	return string

def splitRow(row : bytes) -> dict[str,bytes]:
	cells = [interpret(cell, seps=[]) for cell in row.split(b"\t")]
	rowDict = dict(zip(COLUMNS, cells))
	if len(cells) > len(COLUMNS):
		rowDict["FORMAT"] = cells[len(COLUMNS)]
		rowDict["SAMPLES"] = cells[len(COLUMNS)+1:]
	elif len(cells) < len(COLUMNS):
		return None
	
	return rowDict

def rowFromBytes(row : bytes) -> RowDict:
	rowDict = splitRow(row)
	if rowDict is None:
		LOGGER.warning("Bad row in VCF file.")
		raise ValueError("Bad row in VCF file.")
	return RowDict(rowDict)

class RowDict(dict):
	CHROM : str
	POS : int
	ID : str
	REF : str
	ALT : list
	QUAL : int
	FILTER : str
	INFO : list
	FORMAT : list
	SAMPLES : list

	@property
	def CHROM(self): return self["CHROM"]
	@property
	def POS(self): return self["POS"]
	@property
	def ID(self): return self["ID"]
	@property
	def REF(self): return self["REF"]
	@property
	def ALT(self): return self["ALT"]
	@property
	def QUAL(self): return self["QUAL"]
	@property
	def FILTER(self): return self["FILTER"]
	@property
	def INFO(self): return self["INFO"]
	@property
	def FORMAT(self): return self["FORMAT"]
	@property
	def SAMPLES(self): return self["SAMPLES"]

	@INFO.setter
	def INFO(self, value : str):
		f = lambda p: (p[0], interpret(p[1]))
		self["INFO"] = dict([f(variable.split("=", maxsplit=1)) for variable in value.split(";")])

	@FORMAT.setter
	def FORMAT(self, value : str):
		self["FORMAT"] = value.split(":")

	@SAMPLES.setter
	def SAMPLES(self, samples : str):
		self["SAMPLES"] = [{k:v for k, v in zip(self["FORMAT"], interpret(sample))} for sample in samples]

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.INFO = self["INFO"]
		self.FORMAT = self["FORMAT"]
		self.SAMPLES = self["SAMPLES"]

class VCFIOWrapper:
	file : TextIO

	def __init__(self, filename, mode):
		LOGGER.debug(f"Opening vcf file: open({filename!r}, {mode!r})")
		self.file = open(filename, mode)
	
	def __del__(self):
		self.file.close()

	def close(self):
		self.file.close()

'''Concerns about future-proofing: diploid genomes give more than one result-column per sample. Looks like: 0|1:[...]
where | separates the base call for each chromosomal pair.'''
class CreateVCF(VCFIOWrapper):
	meta : list[str]
	header : list[str] # Header for the VCF, contains names of samples and other columns

	def __init__(self, filename : str, referenceFile, newline : str="\n"):
		super().__init__(filename, mode="w")
		self.newline = newline

		self.header = VCF_HEADER.format(dateYYYYMMDD="{:0>4}{:0>2}{:0>2}".format(*(time.localtime()[:3])), refPath=referenceFile)
		self.file.write(self.header)
			
	def add(self, CHROM : str=".", POS : str=".", ID : str=".", REF : str=".", ALT : str=".", QUAL : str=".", FILTER : str=".", INFO : str="."):
		self.file.write( f"{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\n")
	

class ReadVCF(VCFIOWrapper):
	entryRows : list[int]
	rowsBySelection : dict[str,dict[str|int, set[int]]]
	columns : list[str] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLES"]

	sampleNames : list[str]
	byteToRow : dict[int,int]

	def __init__(self, filename : str):
		super().__init__(filename, mode="rb")

		self.entryRows = []
		self.byteToRow = {}
		self.rowsBySelection = {c:{} for c in self.columns}

		startOfRow = 0
		rowNumber = 0
		# Pass through Header
		for row in self.file:
			rowLength = len(row)
			rowNumber += 1
			self.byteToRow[startOfRow] = rowNumber
			if row.startswith(b"##"):
				pass # Meta row
			elif row.startswith(b"#"):
				l = len("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO")
				if len(row) >= l and len(row) <= l+3:
					self.sampleNames = []
				else:
					self.sampleNames = row[l+len("FORMAT	"):].strip().decode("utf-8").split("\t")
				startOfRow += rowLength
				break

			startOfRow += rowLength

		# Start Indexing the entries.
		for row in self.file:
			rowLength = len(row)
			rowNumber += 1
			self.byteToRow[startOfRow] = rowNumber
			self.entryRows.append(startOfRow)

			rowDict = splitRow(row)

			if rowDict is None:
				LOGGER.error("Bad row in VCF file '{filename}' Row #{n}".format(filename=filename, n=rowNumber))
				raise ValueError("Bad row in VCF file '{filename}' Row #{n}".format(filename=filename, n=rowNumber))
			
			for key in self.columns:
				if key in ["INFO", "FORMAT", "SAMPLES"]: continue

				if key in rowDict:
					try:
						if rowDict[key] in self.rowsBySelection[key]:
							self.rowsBySelection[key][rowDict[key]].add(startOfRow)
						else:
							self.rowsBySelection[key][rowDict[key]] = {startOfRow}
					except TypeError:
						# list is not hashable, skip those for now.
						pass

			startOfRow += rowLength

	def __iter__(self):
		self.file.seek(self.entryRows[0])
		return (rowFromBytes(row.rstrip()) for row in self.file)

	def where(self, CHROM : str | list[str]=None, POS : int|list[int]=None, REF : str | list[str]=None, ALT : str | list[str]=None, QUAL : int | list[int]=None, FILTER : str | list[str]=None, INFO : str | list[str]=None, FORMAT=None) -> list[RowDict]:
		"""Gets dictionary of VCF row values, interpreted into pythonic types as well as it can. Dictionary can be
		queried using keys (or attributes) of the same names as the keyword-arguments for this method."""
		
		flags = {c:locals()[c] for c in self.columns if locals()[c] is not None}
		if len(flags) == 0:
			LOGGER.debug("ReadVCF.where() not given any keyword arguments to search by.")
			return None
		
		# Needs to find intersect of row sets.
		sets : list[set] = []
		for flag, value in flags.items():
			if type(value) in [list, tuple]:
				s = set()
				for key in value:
					s |= self.rowsBySelection[flag][key]
				sets.append(s)
			else:
				sets.append(self.rowsBySelection[flag][value])
		
		try:
			rowStarts = sorted(list(set.intersection(*sets)))
		except:
			LOGGER.warning("ReadVCF.where() was unable to find any rows that match the given query.")
			rowStarts = []
		
		rows : list[RowDict] = []
		for rowStart in rowStarts:
			self.file.seek(rowStart)
			try:
				rows.append(rowFromBytes(self.file.readline().rstrip()))
			except:
				LOGGER.error("Bad row in VCF file '{filename}' Row #{n}".format(filename=self.file.name, n=self.byteToRow[rowStart]))
				raise ValueError("Bad row in VCF file '{filename}' Row #{n}".format(filename=self.file.name, n=self.byteToRow[rowStart]))

		return rows

@overload
def openVCF(filename : str, mode : str, referenceFile : None=None) -> ReadVCF: pass
@overload
def openVCF(filename : str, mode : str, referenceFile : str=None) -> CreateVCF: pass

def openVCF(filename : str, mode : str, referenceFile : str=None) -> ReadVCF|CreateVCF:
	if mode == "r":
		return ReadVCF(filename=filename)
	elif mode == "a":
		raise NotImplementedError("Appending/building on a .vcf file is not yet implemented!")
	elif mode == "w":
		if referenceFile is None:
			raise TypeError("Missing keyword argument 'referenceFile' required to create a .vcf")
		return CreateVCF(filename=filename, referenceFile=referenceFile)
	else:
		raise ValueError(f"openVCF: {mode!r} is not a recognized file mode.")

if __name__ == "__main__":
	vcf = ReadVCF("specCallData.vcf")
	print(f"{vcf.byteToRow=}")
	print(f"{vcf.sampleNames=}")
	print(f"{vcf.rowsBySelection=}")

	for SNP in vcf:
		print(SNP.keys(), *["{}={}".format(key, SNP[key]) for key in COLUMNS[:-1]])