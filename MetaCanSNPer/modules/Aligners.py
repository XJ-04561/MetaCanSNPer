

import logging
import LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)
from Wrappers import Aligner

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''

class ProgressiveMauve(Aligner):
	softwareName = "progressiveMauve"
	commandTemplate = "progressiveMauve {0[options]} {0[ref]} {0[query]} > {0[output]}"
	outFormat = "xmfa"

	'''This is probably not needed.'''
	boolFlags = [
		"--disable-backbone",               # Disable backbone detection
		"--mums",                           # Find MUMs only, do not attempt to determine locally collinear blocks (LCBs)
		"--version",                        # Display software version information
		"--debug",                          # Run in debug mode (perform internal consistency checks–very slow)
		"--collinear",                      # Assume that input sequences are collinear–they have no rearrangements
		"--no-weight-scaling",              # Don’t scale LCB weights by conservation distance and breakpoint distance
		"--skip-refinement",                # Do not perform iterative refinement
		"--skip-gapped-alignment",          # Do not perform gapped alignment
		"--mem-clean",                      # Set this to true when debugging memory allocations
		"--seed-family"                     # Use a family of spaced seeds to improve sensitivity
	]
	valueFlags = [
		"--apply-backbone",                	# <file> # Read an existing sequence alignment in XMFA format and apply backbone statistics to it
		"--seed-weight",                    # <number> # Use the specified seed weight for calculating initial anchors
		"--output",                         # <file> # Output file name. Prints to screen by default
		"--backbone-output",                # <file> # Backbone output file name (optional).
		"--match-input",                    # <file> # Use specified match file instead of searching for matches
		"--input-id-matrix",                # <file> # An identity matrix describing similarity among all pairs of input sequences/alignments
		"--max-gapped-aligner-length",      # <number> # Maximum number of base pairs to attempt aligning with the gapped aligner
		"--input-guide-tree",               # <file> # A phylogenetic guide tree in NEWICK format that describes the order in which sequences will be aligned
		"--output-guide-tree",              # <file> # Write out the guide tree used for alignment to a file
		"--scratch-path-1",                 # <path> # Designate a path that can be used for temporary data storage. Two or more paths should be specified.
		"--scratch-path-2",                 # <path> # Designate a path that can be used for temporary data storage. Two or more paths should be specified.
		"--scoring-scheme",                 # <ancestral|sp_ancestral|sp> # Selects the anchoring score function. Default is extant sum-of-pairs (sp).
		"--max-breakpoint-distance-scale",  # <number [0,1]> # Set the maximum weight scaling by breakpoint distance. Defaults to 0.9
		"--conservation-distance-scale",    # <number [0,1]> # Scale conservation distances by this amount. Defaults to 1
		"--bp-dist-estimate-min-score",     # <number> # Minimum LCB score for estimating pairwise breakpoint distance
		"--gap-open",                       # <number> # Gap open penalty
		"--gap-extend",                     # <number> # Gap extend penalty
		"--substitution-matrix",            # <file> # Nucleotide substitution matrix in NCBI format
		"--weight",                         # <number> # Minimum pairwise LCB score
		"--min-scaled-penalty",             # <number> # Minimum breakpoint penalty after scaling the penalty by expected divergence
		"--hmm-p-go-homologous",            # <number> # Probability of transitioning from the unrelated to the homologous state [0.0001]
		"--hmm-p-go-unrelated"              # <number> # Probability of transitioning from the homologous to the unrelated state [0.000001]
	]

	def handleRetCode(self, returncode : int, prefix : str=""):
		if returncode == 0:
			LOGGER.debug("{prefix}progressiveMauve finished with exitcode 0.".format(prefix=prefix))
		elif returncode == 11:
			LOGGER.warning("{prefix}WARNING progressiveMauve finished with a exitcode: {returncode}\nThis progressiveMauve error is showing up for bad genomes containing short repetitive contigs or sequence contains dashes.".format(prefix=prefix, returncode=returncode))
		elif returncode == -6:
			LOGGER.warning("{prefix}Input sequence is not free of gaps, replace gaps with N and retry!!".format(prefix=prefix))
		else:
			LOGGER.warning("{prefix}WARNING progressiveMauve finished with a non zero exitcode: {returncode}\nThe script will terminate when all processes are finished read {log} for more info".format(prefix=prefix, log=self.logFile,returncode=returncode))

def get(softwareName):
	for c in Aligner.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None