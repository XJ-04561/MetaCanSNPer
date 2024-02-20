
import os
import logging
try:
	import LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	from Wrappers import SNPCaller
	from VCFhandler import OpenVCF, DEFAULT_FORMAT
except:	
	import CanSNPer2.modules.LogKeeper as LogKeeper

	LOGGER = LogKeeper.createLogger(__name__)
	from CanSNPer2.modules.Wrappers import SNPCaller
	from CanSNPer2.modules.VCFhandler import OpenVCF, DEFAULT_FORMAT

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
	commandTemplate = "gatk Mutect2 -R {0[refPath]} -I {0[indexPath]} -L {0[SNPs]} --alleles {0[SNPs]} -O {0[output]} > {0[logFile]}"

	def preProcess(self, data, force : bool=False):
		SNPFiles = {}
		for reference, SNPs in data:
			filename = "{ref}.vcf".format(ref=reference)
			if force is True or not os.path.exists(filename):
				
				f = openVCF(filename, "w", os.path.basename(reference))
				positions = list(SNPs.keys())
				positions.sort()
				for pos in positions:
					p, ref, alt, id = POS=SNPs[pos]
					f.append(CHROM=, POS=p, ID=id, REF=ref, ALT=alt)
				f.save()
			SNPFiles[os.splitext(os.path.basename(reference))[1]] = filename
		self.Lib.setSNPs(SNPFiles)

class GATK_HaplotypeCaller(SNPCaller):
	softwareName = "gatk_HaplotypeCaller"
	commandTemplate = "gatk --java-options '-Xmx4g' HaplotypeCaller -R {ref} -I {target} -alleles {ref}.vcf -O {out}.vcf.gz"

def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None