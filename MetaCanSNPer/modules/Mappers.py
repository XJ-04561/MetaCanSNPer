

import logging
import LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)
from Wrappers import Mapper

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''

class Minimap2(Mapper):
	softwareName = "minimap2"
	commandTemplate = "minimap2 -ax {0[preset]} {0[options]} -R '@RG\tID:{0[queryName]}\tSM:{0[queryName]}' --junc-bed '{0[refPath]}.bed' --sam-hit-only '{0[refPath]}' {0[query]} | samtools sort -T '{0[tmpDir]}' -o {0[output]} && samtools index {0[output]} samtools"
	outFormat = "sam"

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	boolFlags = [
		# Indexing options

		"-H",                # Use homopolymer-compressed (HPC) minimizers. An HPC sequence is constructed by contracting homopolymer runs to a single base. An HPC minimizer is a minimizer on the HPC sequence.
		"--idx-no-seq",      # Don’t store target sequences in the index. It saves disk space and memory but the index generated with this option will not work with -a or -c. When base-level alignment is not requested, this option is automatically applied.

		# Mapping options

		"-D",                # If query sequence name/length are identical to the target name/length, ignore diagonal anchors. This option also reduces DP-based extension along the diagonal.
		"-P",                # Retain all chains and don’t attempt to set primary chains. Options -p and -N have no effect when this option is in use.
		"-X",                # Equivalent to ’-DP --dual=no --no-long-join’. Primarily used for all-vs-all read overlapping.
		"--hard-mask-level", # Honor option -M and disable a heurstic to save unmapped subsequences and disables --mask-len.
		"--no-long-join",    # Disable the long gap patching heuristic. When this option is applied, the maximum alignment gap is mostly controlled by -r.
		"--splice", 	       # Enable the splice alignment mode.
		"--sr", 	           # Enable short-read alignment heuristics. In the short-read mode, minimap2 applies a second round of chaining with a higher minimizer occurrence threshold if no good chain is found. In addition, minimap2 attempts to patch gaps between seeds with ungapped alignment.
		"--for-only", 	     # Only map to the forward strand of the reference sequences. For paired-end reads in the forward-reverse orientation, the first read is mapped to forward strand of the reference and the second read to the reverse stand.
		"--rev-only",        # Only map to the reverse complement strand of the reference sequences.
		"--no-pairing",      # Treat two reads in a pair as independent reads. The mate related fields in SAM are still properly populated.

		# Alignment options

		"--no-end-flt",      # Don’t filter seeds towards the ends of chains before performing base-level alignment.

		# Input/output options

		"-a",                # Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
		"-Q",                # Ignore base quality in the input file.
		"-L",                # Write CIGAR with >65535 operators at the CG tag. Older tools are unable to convert alignments with >65535 CIGAR ops to BAM. This option makes minimap2 SAM compatible with older tools. Newer tools recognizes this tag and reconstruct the real CIGAR in memory.
		"-y",                # Copy input FASTA/Q comments to output.
		"-c",                # Generate CIGAR. In PAF, the CIGAR is written to the ‘cg’ custom tag.
		"--MD",              # Output the MD tag (see the SAM spec).
		"--eqx",             # Output =/X CIGAR operators for sequence match/mismatch.
		"-Y",                # In SAM output, use soft clipping for supplementary alignments.
		"-2",                # Use two I/O threads during mapping. By default, minimap2 uses one I/O thread. When I/O is slow (e.g. piping to gzip, or reading from a slow pipe), the I/O thread may become the bottleneck. Apply this option to use one thread for input and another thread for output, at the cost of increased peak RAM.
		"--paf-no-hit",      # In PAF, output unmapped queries; the strand and the reference name fields are set to ‘*’. Warning: some paftools.js commands may not work with such output for the moment.
		"--sam-hit-only",    # In SAM, don’t output unmapped reads.
		"--version",         # Print version number to stdout

		# Miscellaneous options

		"--no-kalloc",       # Use the libc default allocator instead of the kalloc thread-local allocator. This debugging option is mostly used with Valgrind to detect invalid memory accesses. Minimap2 runs slower with this option, especially in the multi-threading mode.
		"--print-qname",     # Print query names to stderr, mostly to see which query is crashing minimap2.
		"--print-seeds"      # Print seed positions to stderr, for debugging only.
	]

	valueFlags = [
		# Indexing options

		"-k",            # INT # 	Minimizer k-mer length [15]
		"-w",            # INT # 	Minimizer window size [2/3 of k-mer length]. A minimizer is the smallest k-mer in a window of w consecutive k-mers.
		"-I",            # NUM # 	Load at most NUM target bases into RAM for indexing [4G]. If there are more than NUM bases in target.fa, minimap2 needs to read query.fa multiple times to map it against each batch of target sequences. NUM may be ending with k/K/m/M/g/G. NB: mapping quality is incorrect given a multi-part index.
		"-d",            # FILE # 	Save the minimizer index of target.fa to FILE [no dump]. Minimap2 indexing is fast. It can index the human genome in a couple of minutes. If even shorter startup time is desired, use this option to save the index. Indexing options are fixed in the index file. When an index file is provided as the target sequences, options -H, -k, -w, -I will be effectively overridden by the options stored in the index file.
		"--alt",         # FILE # 	List of ALT contigs [null]
		"--alt-drop",    # FLOAT # Drop ALT hits by FLOAT fraction when ranking and computing mapping quality [0.15]

		# Mapping options

		"-f",                # FLOAT|INT1[,INT2] # If fraction, ignore top FLOAT fraction of most frequent minimizers [0.0002]. If integer, ignore minimizers occuring more than INT1 times. INT2 is only effective in the --sr or -xsr mode, which sets the threshold for a second round of seeding.
		"-U",                # INT1[,INT2] # Lower and upper bounds of k-mer occurrences [10,1000000]. The final k-mer occurrence threshold is max{INT1, min{INT2, -f}}. This option prevents excessively small or large -f estimated from the input reference. Available since r1034 and deprecating --min-occ-floor in earlier versions of minimap2.
		"--q-occ-frac",      # FLOAT # Discard a query minimizer if its occurrence is higher than FLOAT fraction of query minimizers and than the reference occurrence threshold [0.01]. Set 0 to disable. Available since r1105.
		"-e",                # INT #	Sample a high-frequency minimizer every INT basepairs [500].
		"-g",                # NUM # 	Stop chain enlongation if there are no minimizers within NUM-bp [10k].
		"-r",                # NUM1[,NUM2] # Bandwidth for chaining and base alignment [500,20k]. NUM1 is used for initial chaining and alignment extension; NUM2 for RMQ-based re-chaining and closing gaps in alignments.
		"-n",                # INT # 	Discard chains consisting of <INT number of minimizers [3]
		"-m",                # INT #	Discard chains with chaining score <INT [40]. Chaining score equals the approximate number of matching bases minus a concave gap penalty. It is computed with dynamic programming.
		"--dual",            # yes|no # If no, skip query-target pairs wherein the query name is lexicographically greater than the target name [yes]
		"-p",                # FLOAT #	Minimal secondary-to-primary score ratio to output secondary mappings [0.8]. Between two chains overlaping over half of the shorter chain (controlled by -M), the chain with a lower score is secondary to the chain with a higher score. If the ratio of the scores is below FLOAT, the secondary chain will not be outputted or extended with DP alignment later. This option has no effect when -X is applied.
		"-N",                # INT #	Output at most INT secondary alignments [5]. This option has no effect when -X is applied.
		"-G",                # NUM #	Maximum gap on the reference (effective with -xsplice/--splice). This option also changes the chaining and alignment band width to NUM. Increasing this option slows down spliced alignment. [200k]
		"-F",                # NUM #	Maximum fragment length (aka insert size; effective with -xsr/--frag=yes) [800]
		"-M",                # FLOAT #	Mark as secondary a chain that overlaps with a better chain by FLOAT or more of the shorter chain [0.5]
		"--rmq",             # no|yes # Use the minigraph chaining algorithm [no]. The minigraph algorithm is better for aligning contigs through long INDELs.
		"--mask-len",        # NUM # Keep an alignment if dropping it leaves an unaligned region on query longer than INT [inf]. Effective without --hard-mask-level.
		"--max-chain-skip",  # INT # A heuristics that stops chaining early [25]. Minimap2 uses dynamic programming for chaining. The time complexity is quadratic in the number of seeds. This option makes minimap2 exits the inner loop if it repeatedly sees seeds already on chains. Set INT to a large number to switch off this heurstics.
		"--max-chain-iter",  # INT # Check up to INT partial chains during chaining [5000]. This is a heuristic to avoid quadratic time complexity in the worst case.
		"--chain-gap-scale", # FLOAT # Scale of gap cost during chaining [1.0]
		"--split-prefix",    # STR # Prefix to create temporary files. Typically used for a multi-part index.
		"--frag",            # no|yes # Whether to enable the fragment mode [no]
		"--heap-sort",       # no|yes # If yes, sort anchors with heap merge, instead of radix sort. Heap merge is faster for short reads, but slower for long reads. [no]

		# Alignment options

		"-A",                # INT #	Matching score [2]
		"-B",                # INT #	Mismatching penalty [4]
		"-O",                # INT1[,INT2] # Gap open penalty [4,24]. If INT2 is not specified, it is set to INT1.
		"-E",                # INT1[,INT2] # Gap extension penalty [2,1]. A gap of length k costs min{O1+k*E1,O2+k*E2}. In the splice mode, the second gap penalties are not used.
		"-C",                # INT #	Cost for a non-canonical GT-AG splicing (effective with --splice) [0]
		"-z",                # INT1[,INT2] # Truncate an alignment if the running alignment score drops too quickly along the diagonal of the DP matrix (diagonal X-drop, or Z-drop) [400,200]. If the drop of score is above INT2, minimap2 will reverse complement the query in the related region and align again to test small inversions. Minimap2 truncates alignment if there is an inversion or the drop of score is greater than INT1. Decrease INT2 to find small inversions at the cost of performance and false positives. Increase INT1 to improves the contiguity of alignment at the cost of poor alignment in the middle.
		"-s",                # INT #	Minimal peak DP alignment score to output [40]. The peak score is computed from the final CIGAR. It is the score of the max scoring segment in the alignment and may be different from the total alignment score.
		"-u",                # CHAR # How to find canonical splicing sites GT-AG - f: transcript strand; b: both strands; n: no attempt to match GT-AG [n]
		"--end-bonus",       # INT # Score bonus when alignment extends to the end of the query sequence [0].
		"--score-N",         # INT # Score of a mismatch involving ambiguous bases [1].
		"--splice-flank",    # yes|no # Assume the next base to a GT donor site tends to be A/G (91% in human and 92% in mouse) and the preceding base to a AG acceptor tends to be C/T [no]. This trend is evolutionarily conservative, all the way to S. cerevisiae (PMID:18688272). Specifying this option generally leads to higher junction accuracy by several percents, so it is applied by default with --splice. However, the SIRV control does not honor this trend (only ~60%). This option reduces accuracy. If you are benchmarking minimap2 on SIRV data, please add --splice-flank=no to the command line.
		"--junc-bed",        # FILE # Gene annotations in the BED12 format (aka 12-column BED), or intron positions in 5-column BED. With this option, minimap2 prefers splicing in annotations. BED12 file can be converted from GTF/GFF3 with ‘paftools.js gff2bed anno.gtf’ [].
		"--junc-bonus",      # INT # Score bonus for a splice donor or acceptor found in annotation (effective with --junc-bed) [9].
		"--end-seed-pen",    # INT # Drop a terminal anchor if s<log(g)+INT, where s is the local alignment score around the anchor and g the length of the terminal gap in the chain. This option is only effective with --splice. It helps to avoid tiny terminal exons. [6]
		"--cap-sw-mem",      # NUM # Skip alignment if the DP matrix size is above NUM. Set 0 to disable [100m].
		"--cap-kalloc",      # NUM # Free thread-local kalloc memory reservoir if after the alignment the size of the reservoir above NUM. Set 0 to disable [0].

		# Input/output options

		"-o",                # FILE # Output alignments to FILE [stdout].
		"-R",                # STR # SAM read group line in a format like @RG\\tID:foo\\tSM:bar [].
		"--cs",              # [=STR] # Output the cs tag. STR can be either short or long. If no STR is given, short is assumed. [none]
		"--seed",            # INT # Integer seed for randomizing equally best hits. Minimap2 hashes INT and read name when choosing between equally best hits. [11]
		"-t",                # INT # Number of threads [3]. Minimap2 uses at most three threads when indexing target sequences, and uses up to INT+1 threads when mapping (the extra thread is for I/O, which is frequently idle and takes little CPU time).
		"-K",                # NUM # Number of bases loaded into memory to process in a mini-batch [500M]. Similar to option -I, K/M/G/k/m/g suffix is accepted. A large NUM helps load balancing in the multi-threading mode, at the cost of increased memory.
		"--secondary",       # yes|no # Whether to output secondary alignments [yes]
		"--max-qlen",        # NUM # Filter out query sequences longer than NUM.

		# Preset options

		"-x"                 # STR # Preset []. This option applies multiple options at the same time. It should be applied before other options because options applied later will overwrite the values set by -x. Available STR are:
			# map-ont   Align noisy long reads of ~10% error rate to a reference genome. This is the default mode.
			# map-hifi 	Align PacBio high-fidelity (HiFi) reads to a reference genome (-k19 -w19 -U50,500 -g10k -A1 -B4 -O6,26 -E2,1 -s200).
			# map-pb    Align older PacBio continuous long (CLR) reads to a reference genome (-Hk19).
			# asm5 	    Long assembly to reference mapping (-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%.
			# asm10 	  Long assembly to reference mapping (-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B9 -O16,41 -E2,1 -s200 -z200 -N50). Up to 10% sequence divergence.
			# asm20 	  Long assembly to reference mapping (-k19 -w10 -U50,500 --rmq -r100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50). Up to 20% sequence divergence.
			# splice 	  Long-read spliced alignment (-k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -b0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes). In the splice mode, 1) long deletions are taken as introns and represented as the ‘N’ CIGAR operator; 2) long insertions are disabled; 3) deletion and insertion gap costs are different during chaining; 4) the computation of the ‘ms’ tag ignores introns to demote hits to pseudogenes.
			# splice:hq Long-read splice alignment for PacBio CCS reads (-xsplice -C5 -O6,24 -B4).
			# sr        Short single-end reads without splicing (-k21 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -b0 -r100 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g100 -2K50m --heap-sort=yes --secondary=no).
			# ava-pb    PacBio CLR all-vs-all overlap mapping (-Hk19 -Xw5 -e0 -m100).
			# ava-ont   Oxford Nanopore all-vs-all overlap mapping (-k15 -Xw5 -e0 -m100 -r2k).
	]

	def preProcess(self, references, snps):
		# Create bed-file of only the SNP locations + padding
		pad = self.settings["RegionPadding"]
		for genome, (refPath,*_) in self.Lib.getReferences().items():
			LOGGER.debug("open('{}'+'.bed', 'w')").format(refPath)
			f = open(refPath+".bed", "w")
			# Write a header maybe??
			for pos in sorted(self.Lib.SNPs):
				f.write("{}\t{}\t{}\n".format(start=pos-pad-2, end=pos+pad))
			f.close()
		if "Preset" in self.settings:
			self.formatDict["preset"] = self.settings["Preset"]
			'''Valid presets:
			map-ont		Align noisy long reads of ~10% error rate to a reference genome. This is the default mode.
			map-hifi	Align PacBio high-fidelity (HiFi) reads to a reference genome (-k19 -w19 -U50,500 -g10k -A1 -B4 -O6,26 -E2,1 -s200).
			map-pb		Align older PacBio continuous long (CLR) reads to a reference genome (-Hk19).
			asm5		Long assembly to reference mapping (-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B19 -O39,81 -E3,1 -s200 -z200 -N50). Typically, the alignment will not extend to regions with 5% or higher sequence divergence. Only use this preset if the average divergence is far below 5%.
			asm10		Long assembly to reference mapping (-k19 -w19 -U50,500 --rmq -r100k -g10k -A1 -B9 -O16,41 -E2,1 -s200 -z200 -N50). Up to 10% sequence divergence.
			asm20		Long assembly to reference mapping (-k19 -w10 -U50,500 --rmq -r100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50). Up to 20% sequence divergence.
			splice		Long-read spliced alignment (-k15 -w5 --splice -g2k -G200k -A1 -B2 -O2,32 -E1,0 -b0 -C9 -z200 -ub --junc-bonus=9 --cap-sw-mem=0 --splice-flank=yes). In the splice mode, 1) long deletions are taken as introns and represented as the ‘N’ CIGAR operator; 2) long insertions are disabled; 3) deletion and insertion gap costs are different during chaining; 4) the computation of the ‘ms’ tag ignores introns to demote hits to pseudogenes.
			splice:hq	Long-read splice alignment for PacBio CCS reads (-xsplice -C5 -O6,24 -B4).
			sr			Short single-end reads without splicing (-k21 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -b0 -r100 -p.5 -N20 -f1000,5000 -n2 -m20 -s40 -g100 -2K50m --heap-sort=yes --secondary=no).
			ava-pb		PacBio CLR all-vs-all overlap mapping (-Hk19 -Xw5 -e0 -m100).
			ava-ont		Oxford Nanopore all-vs-all overlap mapping (-k15 -Xw5 -e0 -m100 -r2k). 
			'''
		else:
			LOGGER.error("No preset for minimap2 set!")
			raise ValueError("No preset for minimap2 set!")
		

def get(softwareName) -> Mapper:
	for c in Mapper.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None