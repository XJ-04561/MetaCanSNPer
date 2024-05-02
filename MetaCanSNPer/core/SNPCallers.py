
from MetaCanSNPer.Globals import *
import MetaCanSNPer.core.LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)

from MetaCanSNPer.core.Wrappers import SNPCaller

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''
class ParseXMFA2(SNPCaller):
	softwareName = "ParseXMFA2"
	commandTemplate = "ParseXMFA2 {alignmentPath!r} {targetSNPs!r} -rID 1 -o {output!r}"
	inFormat = ["xmfa"]
	outFormat = "vcf"

class GATK_Mutect2(SNPCaller):
	softwareName = "gatk_Mutect2"
	commandTemplate = "samtools faidx {refPath!r} && gatk CreateSequenceDictionary -R {refPath!r}; gatk IndexFeatureFile -I {targetSNPs!r} --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' && gatk Mutect2 --genotype-germline-sites --genotype-pon-sites --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -R {refPath!r} -I {mapPath!r} -L {targetSNPs!r} --force-call-filtered-alleles --alleles {targetSNPs!r} -O {output!r}"
	inFormat = ["bam"]
	outFormat = "vcf"

class GATK_HaplotypeCaller(SNPCaller):
	softwareName = "gatk_HaplotypeCaller"
	commandTemplate = "gatk IndexFeatureFile -I {targetSNPs!r} --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' && gatk HaplotypeCaller --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' -R {refPath!r} -I {mapPath!r} -L {targetSNPs!r} --alleles {targetSNPs!r} -O {output!r}"
	inFormat = ["bam"]
	outFormat = "vcf"

def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName.lower() == softwareName.lower():
			return c
	return None