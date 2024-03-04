
import os
import time
try:
	import LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	from Wrappers import SNPCaller
	from Databases import DatabaseReader
except:
	import MetaCanSNPer.modules.LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	from MetaCanSNPer.modules.Wrappers import SNPCaller
	from MetaCanSNPer.modules.Databases import DatabaseReader

VCF_HEADER = """##fileformat=VCFv4.3
##fileDate={dateYYYYMMDD}
##source=MetaCanSNPer
##reference={refPath}
##contig=<ID=0,URL={refPath}>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"""
VCF_ROW =       "{CHROM}	{POS}	{ID}	{REF}	{ALT}	{QUAL}	{FILTER}	{INFO}"

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''

class GATK_Mutect2(SNPCaller):
	softwareName = "gatk_Mutect2"
	commandTemplate = "gatk IndexFeatureFile -I '{0[SNPs]}' && gatk Mutect2 -R {0[refPath]} -I {0[indexPath]} -L {0[SNPs]} --alleles {0[SNPs]} -O {0[output]} > {0[logFile]}"

class GATK_HaplotypeCaller(SNPCaller):
	softwareName = "gatk_HaplotypeCaller"
	commandTemplate = "gatk IndexFeatureFile -I '{0[SNPs]}' && gatk HaplotypeCaller -R {0[refPath]} -I {0[indexPath]} -L {0[SNPs]} --alleles {0[SNPs]} -O {0[output]} > {0[logFile]}"

def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None