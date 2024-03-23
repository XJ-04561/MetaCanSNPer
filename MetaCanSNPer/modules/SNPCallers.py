
from MetaCanSNPer.Globals import *
import MetaCanSNPer.modules.LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)

from MetaCanSNPer.modules.Wrappers import SNPCaller

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''
class ParseXMFA2(SNPCaller):
	softwareName = "ParseXMFA2"
	commandTemplate = "ParseXMFA2 {alignmentPath!r} {targetSNPs!r} -rID 1 -o {output!r} > {logFile!r}"
	inFormat = ["xmfa"]
	outFormat = "vcf"

class GATK_Mutect2(SNPCaller):
	softwareName = "gatk_Mutect2"
	commandTemplate = "gatk IndexFeatureFile -I {targetSNPs!r} && gatk Mutect2 --genotype-germline-sites --genotype-pon-sites -R {refPath!r} -I {mapPath!r} -L {targetSNPs!r} --force-call-filtered-alleles --alleles {targetSNPs!r} -O {output!r} > {logFile!r}"
	inFormat = ["bam"]
	outFormat = "vcf"

class GATK_HaplotypeCaller(SNPCaller):
	softwareName = "gatk_HaplotypeCaller"
	commandTemplate = "gatk IndexFeatureFile -I {targetSNPs!r} && gatk HaplotypeCaller -R {refPath!r} -I {mapPath!r} -L {targetSNPs!r} --alleles {targetSNPs!r} -O {output!r} > {logFile!r}"
	inFormat = ["bam"]
	outFormat = "vcf"

def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName.lower() == softwareName.lower():
			return c
	return None