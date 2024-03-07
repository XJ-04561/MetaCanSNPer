
import sys
try:
	import MetaCanSNPer.modules.LogKeeper as LogKeeper
	from MetaCanSNPer.modules.Wrappers import SNPCaller
except:
	import LogKeeper as LogKeeper
	from Wrappers import SNPCaller

	LOGGER = LogKeeper.createLogger(__name__)
PYTHON_INTERPRETER = sys.executable

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''
class ParseXMFA2(SNPCaller):
	softwareName = "ParseXMFA2"
	commandTemplate = "ParseXMFA2 '{0[alignmentPath]}' '{0[targetSNPs]}' -rID 1 -o {0[output]} > {0[logFile]}"
	inFormat = [".xmfa"]
	outFormat = ".vcf"

class GATK_Mutect2(SNPCaller):
	softwareName = "gatk_Mutect2"
	commandTemplate = "gatk IndexFeatureFile -I '{0[targetSNPs]}' && gatk Mutect2 -R {0[refPath]} -I {0[mapPath]} -L {0[targetSNPs]} --force-call-filtered-alleles --alleles {0[targetSNPs]} -O {0[output]} > {0[logFile]}"
	inFormat = [".bam"]
	outFormat = ".vcf"

class GATK_HaplotypeCaller(SNPCaller):
	softwareName = "gatk_HaplotypeCaller"
	commandTemplate = "gatk IndexFeatureFile -I '{0[targetSNPs]}' && gatk HaplotypeCaller -R {0[refPath]} -I {0[mapPath]} -L {0[targetSNPs]} --alleles {0[targetSNPs]} -O {0[output]} > {0[logFile]}"
	inFormat = [".bam"]
	outFormat = ".vcf"

def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None