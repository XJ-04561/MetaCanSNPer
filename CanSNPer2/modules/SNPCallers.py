

import logging
import LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)
from Wrappers import SNPCaller

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''

class Timeout(SNPCaller):
	softwareName = "timeout"
	commandTemplate = "timeout {ref} {target} > {output}"

class Sleep(SNPCaller):
	softwareName = "sleep"
	commandTemplate = "sleep {ref} {target} > {output}"

class GATK_Mutect2(SNPCaller):
	softwareName = "gatk_Mutect2"
	commandTemplate = "gatk Mutect2 -R {ref} -I {target} --alleles {ref}.vcf -O {out}.vcf.gz"

class GATK_HaplotypeCaller(SNPCaller):
	softwareName = "gatk_HaplotypeCaller"
	commandTemplate = "gatk --java-options '-Xmx4g' HaplotypeCaller -R {ref} -I {target} -alleles {ref}.vcf -O {out}.vcf.gz"

def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None