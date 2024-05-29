'''
	All that is needed to create a new implementation is to inherit from the correct software type (`Aligner` in this case) and set
	the two class attributes accordingly.
'''

from MetaCanSNPer.Globals import *

from MetaCanSNPer.core.Wrappers import Mapper

class Minimap2(Mapper):
	softwareName = "minimap2"
	commandTemplate = "minimap2 -a {options} -R '@RG\\tID:{queryName}\\tSM:{queryName}' --sam-hit-only {refPath!r} {query} -o '{output}.sam' && samtools sort --reference {refPath!r} -O BAM -T {refName!r} -o {output!r} '{output}.sam' && samtools index {output!r} && rm -f '{output}.sam'"
	outFormat = "bam"

	boolFlags = [
		"-H", "--idx-no-seq", "-D", "-P", "-X", "--hard-mask-level", "--no-long-join", "--splice", "--sr",
		"--for-only", "--rev-only", "--no-pairing", "--no-end-flt", "-a", "-Q", "-L", "-y", "-c", "--MD", "--eqx",
		"-Y", "-2", "--paf-no-hit", "--sam-hit-only", "--version", "--no-kalloc", "--print-qname", "--print-seeds"
	]

	valueFlags = [
		"-k", "-w", "-I", "-d", "--alt", "--alt-drop", "-f", "-U", "--q-occ-frac", "-e", "-g", "-r", "-n", "-m",
		"--dual", "-p", "-N", "-G", "-F", "-M", "--rmq", "--mask-len", "--max-chain-skip", "--max-chain-iter",
		"--chain-gap-scale", "--split-prefix", "--frag", "--heap-sort", 	"-A", "-B", "-O", "-E", "-C", "-z",
		"-s", "-u", "--end-bonus", "--score-N", "--splice-flank", "--junc-bed", "--junc-bonus", "--end-seed-pen",
		"--cap-sw-mem", "--cap-kalloc", "-o", "-R", "--cs", "--seed", "-t", "-K", "--secondary", "--max-qlen", "-x"
	]

	def preProcess(self):
		pass
		

def get(softwareName) -> Mapper:
	for c in Mapper.__subclasses__():
		if c.softwareName.lower() == softwareName.lower():
			return c
	return None