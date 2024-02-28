
import numpy as np
import time

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

VCF_ROW =       "{CHROM}	{POS}	{ID}	{REF}	{ALT}	{QUAL}	{FILTER}	{INFO}"

'''Concerns about future-proofing: diploid genomes give more than one result-column per sample. Looks like: 0|1:[...]
where | separates the base call for each chromosomal pair.'''
class CreateVCF:
	meta : list[str]
	header : list[str] # Header for the VCF, contains names of samples and other columns
	positions : np.ndarray[int] # List of positions of the SNPs
	SNPsinfo : np.ndarray[str] # A lookup for information about the SNP at that position
	samples : np.ndarray[np.ndarray] # A 3D table of the data. axis 0 separates positions, axis 1 separates samples, axis 2 separates data
	columnIndex : dict[str, int]= {
		"CHROM" : 0,
		"POS" : 1,
		"ID" : 2,
		"REF" : 3,
		"ALT" : 4,
		"QUAL" : 5,
		"FILTER" : 6,
		"INFO" : 7,
		"FORMAT" : 8
	}

	def __init__(self, filename : str, referenceFile, newline : str="\n"):
		self.fileHandler = open(filename, mode="w")
		self.newline = newline

		self.header = VCF_HEADER.format(dateYYYYMMDD="{:0>4}{:0>2}{:0>2}".format(*(time.localtime()[:3])), refPath=referenceFile)
		self.fileHandler.write(self.header)
			
	def append(self, CHROM : str=".", POS : str=".", ID : str=".", REF : str=".", ALT : str=".", QUAL : str=".", FILTER : str=".", INFO : str=".", FORMAT : str="."):
		self.fileHandler.write( VCF_ROW.format(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)+self.newline)
	
	def close(self):
		self.fileHandler.close()