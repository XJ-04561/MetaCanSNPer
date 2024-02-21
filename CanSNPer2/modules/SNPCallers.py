
import os
import time
try:
	import LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	from Wrappers import SNPCaller
except:	
	import CanSNPer2.modules.LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	from CanSNPer2.modules.Wrappers import SNPCaller

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

class Timeout(SNPCaller):
	softwareName = "timeout"
	commandTemplate = "timeout {0[options]} {0[ref]} {0[query]} > {0[output]}"

class Sleep(SNPCaller):
	softwareName = "sleep"
	commandTemplate = "sleep {0[options]} {0[ref]} {0[query]} > {0[output]}"

class GATK_Mutect2(SNPCaller):
	softwareName = "gatk_Mutect2"
	commandTemplate = "gatk IndexFeatureFile -I '{0[SNPs]}' && gatk Mutect2 -R {0[refPath]} -I {0[indexPath]} -L {0[SNPs]} --alleles {0[SNPs]} -O {0[output]} > {0[logFile]}"

	def preProcess(self, data, force : bool=False):
		
		# Create VCF files that contain the to-be called SNPs

		SNPFiles = {}
		for genome, SNPs in data:
			refPath, strain, genkbank_id, refseq_id, assembly_name = self.Lib.references[genome] # dadsadasd
			# accession = open(refPath, "r").readline()[1:].split()[0]
			filename = "{ref}.vcf".format(ref=refPath)
			if force is True or not os.path.exists(filename):
				
				vcfFile = open(filename, "w")
				vcfFile.write(VCF_HEADER.format(dateYYYYMMDD="{:0>4}{:0>2}{:0>2}".format(*(time.localtime()[:3])), refPath=refPath))
				positions = list(SNPs.keys())
				positions.sort()
				for pos in positions:
					p, ref, alt, id =SNPs[pos]
					# CHROM has to be the same as the accession id that is in the reference file.
					vcfFile.write(VCF_ROW.format(CHROM="0", POS=p, ID=id, REF=ref, ALT=alt)+"\n")
				vcfFile.close()
				self.commandTemplate = "gatk IndexFeatureFile -I '{0[SNPs]}' && " + self.commandTemplate
				
			SNPFiles[genome] = filename
		self.Lib.setSNPfiles(SNPFiles)

class GATK_HaplotypeCaller(SNPCaller):
	softwareName = "gatk_HaplotypeCaller"
	commandTemplate = "gatk --java-options '-Xmx4g' HaplotypeCaller -R {ref} -I {target} -alleles {ref}.vcf -O {out}.vcf.gz"

def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None