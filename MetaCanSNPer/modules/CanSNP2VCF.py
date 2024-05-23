import os

def CanSNP2VCF(lib : "DirectoryLibrary", force : bool=False):
	from MetaCanSNPer.modules.Database import Chromosome, Position, GenomeID
	from VariantCallFixer import openVCF
	for genomeID, genome, refPath, SNPPath in map(lambda x: (x[0], x[1], lib.references[x[1]], lib.targetSNPs[x[1]]), lib.database.references):
		if SNPPath is not None and not force:
			continue
		
		filename = f"{genome}.vcf"
		path = lib.SNPDir.writable

		if path is None:
			raise PermissionError("Path(s) are not writable, so no `.VCF`-file can be created.\n" f"Path: {lib.SNPDir!r}")
		if not os.path.exists(refPath):
			raise FileNotFoundError(f"Could not find file for reference genome {genome!r}. Either the file has not been downloaded or the path to the file is not correct. Path to file was {refPath!r}")

		filePath = path / filename
		
		with openVCF(filePath, mode="w", referenceFile=refPath) as vcfFile:
			for chromosome, pos in lib.database[Chromosome, Position, GenomeID == genomeID]:
				vcfFile.add(CHROM=chromosome, POS=pos, REF="N", ALT="A,T,C,G")
			lib.targetSNPs[genome] = filePath

try:
	from MetaCanSNPer.core.DirectoryLibrary import DirectoryLibrary
except:
	pass