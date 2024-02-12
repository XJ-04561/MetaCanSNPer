'''
CanSNPer2 module: A toolkit for SNP-typing using NGS data.
Copyright (C) 2019 David Sundell @ FOI bioinformatics group
'''

import os
import logging
try:
	from CanSNPer2.modules.LogKeeper import createLogger

	LOGGER = createLogger(__name__)

	## import CanSNPer2 specific modules
	from CanSNPer2.modules.ParseXMFA2 import ParseXMFA
	from CanSNPer2.modules.DatabaseConnection import CanSNPdbFunctions
	from CanSNPer2.modules.NewickTree import NewickTree
	from CanSNPer2.CanSNPerTree import __version__
	from CanSNPer2.modules.DirectoryLibrary import DirectoryLibrary
	from CanSNPer2.modules.Wrappers import Aligner, Mapper, SNPCaller
	import CanSNPer2.modules.Aligners as Aligners
	import CanSNPer2.modules.Mappers as Mappers
	import CanSNPer2.modules.SNPCallers as SNPCallers
except:
	from LogKeeper import createLogger

	LOGGER = createLogger(__name__)

	## import CanSNPer2 specific modules
	from ParseXMFA2 import ParseXMFA
	from DatabaseConnection import CanSNPdbFunctions
	from NewickTree import NewickTree
	from DirectoryLibrary import DirectoryLibrary
	from Wrappers import Aligner, Mapper, SNPCaller
	import Aligners
	import Mappers
	import SNPCallers

## import standard python libraries for subprocess and multiprocess
from subprocess import Popen,PIPE,STDOUT
from multiprocessing import Process, Queue
from time import sleep,time
try:
	import tomllib as toml
except ModuleNotFoundError:
	try:
		import tomli as toml
	except:
		raise ModuleNotFoundError("TOML-reading module not found. For Python 3.5> it should be included as 'tomllib', if older than 2.5, then install 'tomli'.")

import random
random.seed()

class Error(Exception):
	"""docstring for Error"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


class CanSNPer2Error(Error):
	"""docstring for MauveE"""
	pass

class CanSNPer2:
	outputTemplate = "{ref}_{target}.{format}"
	databaseName : str
	database : CanSNPdbFunctions
	Lib : DirectoryLibrary[str]
	settings : dict

	"""docstring for CanSNPer2"""
	def __init__(self, tmpDir=None, refDir="References", database="CanSNPer.fdb", settingsFile : str=None, **kwargs):
		
		self.databaseName = database
		self.Lib = DirectoryLibrary(tmpDir=tmpDir, refDir=refDir)
		settings : dict = toml.load(open(self.Lib.get("defaultFlags.toml", "installDir") if settingsFile is None else self.Lib.get(settingsFile, "workDir")), "rb")
		# Settings hierarchy looks like this: ["Category"]["Flag"] -> Value

		# Flatten hierarchy so flags are easily accessible
		self.settings = {}
		for flags in settings.values():
			for flag, value in flags.items():
				if flag in kwargs:
					self.settings[flag] = kwargs[flag]
				else:
					self.settings[flag] = value

		# Need for this is unknown, might remove later.
		self.summarySet = set()
		self.calledGenome = {}

	'''Database functions'''

	def setDatabase(self, database : str):
		path = self.Lib.get(database, "databaseDir")
		self.databaseName = path

	def connectDatabase(self):
		self.database = CanSNPdbFunctions(self.databaseName)

	'''MetaCanSNPer set functions'''

	def setQuery(self, query : str):
		self.Lib.setQuery(query)
	
	def setRefDir(self, references : str):
		''''references' is either an absolut path or a path to an accepted directory.'''
		self.Lib.setRefDir(references)
	
	def setOutDir(self, outDir : str):
		self.Lib.setOutDir(outDir)

	'''MetaCanSNPer get functions'''

	def getTmps(self):
		'''Returns a dictionary of all files in and under the temporary directory, where the keys are the extensions of
		the files, and the values are lists of the paths to all files with that extension.'''
		files = []
		for dirpath, dirnames, filenames in os.walk(self.Lib.tmpDir):
			files.extend([os.path.join(dirpath, f) for f in filenames])
		formats = [f.rsplit(".")[-1] for f in files]
		return {format:[f for f in files if f.endswith(format)] for format in formats}

	def getReferences(self):
		localReferences = self.Lib.getReferences()
		localReferenceIDs = map(lambda f : os.path.basename(f).rsplit(".")[0], localReferences)
		for ref, refID in self.database.get_genomes():
			# refID should be the genbank ID of the given genome.
			if refID in localReferenceIDs:
				# genome exists in the filesystem.
			else:

		## DownloadGenomes - THIS DOES WHAT YOU NEED

	def getQuery(self):
		return self.Lib.query

	'''Indexing methods'''

	def createIndex(self, software : str=None, kwargs : dict={}):
		'''Align sequences using subprocesses.'''

		LOGGER.debug("Checking for Aligner named '{}'".format( software))
		indexer : Aligner = Aligners.get(software) # Fetch the object class of the specified aligner software.
		if indexer is None:
			LOGGER.debug("Checking for Mapper named '{}'".format( software))
			indexer : Mapper = Mappers.get(software) # Fetch the object class of the specified Mapper software.
		else:
			# No Aligner or Mapper found that is implemented yet.
			LOGGER.error("No software defined for name '{}'".format( software))
			raise NotImplementedError("No software defined for name '{}'".format( software))
		
		indexer(self.Lib, self.outputTemplate, kwargs=kwargs)

		output = indexer.start()

		# Start aligning
		while not indexer.finished():
			finished = indexer.waitNext()
			for i in finished:
				key, path = output[i][0]
				if key not in self.Lib.indexed and indexer.returncodes[i][-1] == 0:
					self.Lib.indexed[key] = path

		# Check that error did not occur.
		while indexer.hickups():
			if not indexer.fixable():
				# Error does not have a known solution that has worked during this run.
				return []
			else:
				# Error has a possible fix, implement the fix.
				indexer.planB()
			
			# Run the alignment just like before, but hopefully fixed.
			while not indexer.finished():
				finished = indexer.waitNext()
				for i in finished:
					key, path = output[i][0]
					if indexer.returncodes[i][-1] == 0:
						self.Lib.indexed[key] = path

	def callSNPs(self, snps : dict[str]):

		for ref, refID in self.database.get_genomes():
			SNPs, snpList = self.database.get_snps(reference=ref)



	'''Functions'''

	def create_tree(self,SNPS,name,called_snps,save_tree,min_required_hits,strictness=0.7, summary=False):
		'''This function uses ETE3 to color the SNP tree in the database with SNPS found in the reference database
			and outputs a pdf file
		'''
		newickTree = NewickTree(self.database,name,self.Lib.outDir,min_required_hits=min_required_hits, strictness=strictness)
		final_snp = newickTree.draw_ete3_tree(SNPS,called_snps,save_tree,summary=summary)
		LOGGER.info("{outdir}/{name}_tree.pdf".format(outdir =self.outdir, name=name))
		return final_snp

	def readQueriesFrom(self, queryFile : str):
		'''If query input is a text file, parse file'''
		with open(queryFile, "r") as f:
			query = [query.strip() for query in f if len(query.strip())]
		return query

	def cleanup(self):
		'''Remove files in temporary folder'''
		LOGGER.info("Clean up temporary files... ")
		del self.Lib
		LOGGER.info("Done!")
		return

	# UNUSED - DELETE?
	def parse_xmfa(XMFA_obj, xmfa_file,results=[]):
		'''Process xmfa file using ParseXMFA object'''
		XMFA_obj.run(xmfa_file)
		results.put(XMFA_obj.get_snps())
		export_results.put(XMFA_obj.get_snp_info())
		called_snps.put(XMFA_obj.get_called_snps())
		return results
	
	def find_snps(self,XMFA_obj,xmfa_file,results=[],export_results=[],called_snps=[]):
		'''Align sequences to references and return SNPs'''
		# execute snp calling
		XMFA_obj.run(xmfa_file)
		#/
		# stream output into Queues (if the complete job output is put in then the queue gets filled when having "large" outputs. It has to be streamed to the queue so that the queue can simultanously be cleared)
		for result in XMFA_obj.get_snps().items(): # returns dictionary
			results.put(result)
		for result in XMFA_obj.get_snp_info(): # returns array
			export_results.put(result)
		for result in XMFA_obj.get_called_snps(): # returns array
			called_snps.put(result)
		#/
		# when finished, put a "stop-signal"
		results.put('XMFA_obj_finish_worker_'+XMFA_obj.reference)
		export_results.put('XMFA_obj_finish_worker_'+XMFA_obj.reference)
		called_snps.put('XMFA_obj_finish_worker_'+XMFA_obj.reference)
		#/

	def find_snps_multiproc(self,xmfa_obj,xmfa_files,export=False):
		'''function to run genomes in paralell'''
		jobs = []
		SNPS = {}
		SNP_info = []
		called_snps = []
		result_queue = Queue()
		export_queue = Queue()
		called_queue = Queue()
		for enum,xmfa_file in enumerate(xmfa_files):
			p = Process(target=self.find_snps, args=(xmfa_obj,xmfa_file ,result_queue,export_queue,called_queue))
			p.start()
			jobs.append(p)
			sleep(0.05) ## A short sleep to make sure all threads do not initiate access to the database file simultanously
		
		## Parse output queues from jobs continuously
		finish_signals = {'SNPS':set(),'SNP_info':set(),'called_snps':set()} # will keep track when all workers are finished
		tic = time() # "tic-toc" clock-sound
		while True:
			# check SNP queue
			if not result_queue.empty():
				tmp_data = result_queue.get()
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['SNPS'].add(tmp_data.split('_')[-1])
					LOGGER.info("Output queue for 'SNPS' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					SNPS[tmp_data[0]] = tmp_data[1]
			#/
			# check SNP info queue
			if not export_queue.empty():
				tmp_data = export_queue.get()
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['SNP_info'].add(tmp_data.split('_')[-1])
					LOGGER.info("Output queue for 'SNP_info' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					SNP_info.append(tmp_data)
			#/
			# check called snps queue
			if not called_queue.empty():
				tmp_data = called_queue.get()
				if type(tmp_data) == str and tmp_data.find('XMFA_obj_finish_worker_') != -1:
					finish_signals['called_snps'].add(tmp_data.split('_')[-1])
					LOGGER.info("Output queue for 'called_snps' finished for reference {reference}".format(reference=tmp_data.split('_')[-1]))
				else:
					called_snps.append(tmp_data)
			#/
			# check if finished (all xmfa files has returned the "finish worker" in all output queues)
			finish_signals_called = []
			for queue_name,signals in finish_signals.items():
				if len(signals) == len(xmfa_files):
					finish_signals_called.append(queue_name)
			if len(finish_signals) == len(finish_signals_called):
				LOGGER.info("All output queues finished!")
				break
			#/
			# check if we have computed for X amount of time, then break and assume something is wrong
			max_time_seconds = 5*60 # 5*60 => 5 minutes
			toc = time()
			time_spent = toc - tic
			if time_spent > max_time_seconds:
				LOGGER.warning("Maximum time spent reached to capture queue output. This is not expected to happen and might mean that your output data is corrupted.")
				break
			#/
		##/
		## Wait until processes terminated (their output should already have been processed in the queue's)
		for job in jobs:
			job.join()
		##/
		return SNPS,SNP_info,called_snps

	def read_result_dir(self):
		'''Fetch all final SNPs from result file after run to produce summary'''
		snplist = set()
		for root, dirs, files in os.walk(self.outdir, topdown=False):
			for f in files:
				if f.endswith("_snps.txt"):
					with open(os.path.join(root,f)) as snpfile:
						for row in snpfile:
							if row.startswith("SNP path:"):
								snp = row.split(";")[-1].strip()
								snplist.add(snp)
								self.called_genome[snp] = f.rsplit(".fasta",1)[0]
		return snplist

	def print_summary(self):
		'''Create summary file and tree'''
		summarypath = "{outdir}/snps_summary.txt".format(outdir=self.outdir)
		SNPS = self.read_result_dir()
		with open(summarypath,"w") as summaryout:
			for snp in SNPS:
				print("{query}: {SNP}".format(query=self.called_genome[snp], SNP=snp),file=summaryout)

		'''Print summary tree showing all unique SNPs in the final tree'''
		self.create_tree([],"summary",SNPS,True,min_required_hits=self.min_required_hits,summary=True)

	def run(self,database):
		'''Run CanSNPer2'''
		LOGGER.info("Running CanSNPer2 version-{version}".format(version=__version__))
		if len(self.query) > 0:
			'''Read query input file if a txt file is supplied insead of fasta files'''
			if self.query[0].endswith(".txt"):
				LOGGER.info("Textfile input was found, parsing filepaths in {q} file".format(q=self.query[0]))
				self.query=self.read_query_textfile_input(self.query)

			'''Main function of CanSNPer2
					1. Align sequences with progressiveMauve
					2. Parse XMFA files and find SNPs
					3. Create a tree visualising SNPs found in sequence
					3b. If requested create a list with SNPs and their status
					4. Clean up tmp directory
					'''

			''' Create ParseXMFA object'''
			parse_xmfa_obj = ParseXMFA(
						database=database,
						export=self.export,
						verbose=self.verbose)  ## Create XMFA object and connect to database
			'''Walk through the list of queries supplied'''
			if not self.skip_mauve: print("Run {n} alignments to references using progressiveMauve".format(n=len(self.query)))
			for q in self.query:			## For each query file_path
				try:
					self.query_name = os.path.basename(q).rsplit(".",1)[0]  ## get name of file and remove ending

					qfile = q.rsplit("/")[-1]   ## Remove path from query name
					if not os.path.exists(q):
						raise FileNotFoundError("Input file: {qfile} was not found!".format(qfile=q))

					outputfile = "{outdir}/{xmfa}_snps.txt".format(outdir=self.outdir,xmfa=self.query_name)
					if os.path.exists(outputfile) and not self.rerun:
						LOGGER.debug("{outputfile} already exits, skip!".format(outputfile=outputfile))
						continue
					LOGGER.info("Running CanSNPer2 on {query}".format(query=qfile))
					if not self.skip_mauve: ### If mauve command was already run before skip step
						LOGGER.info("Run mauve alignments")

					'''For each query fasta align to all CanSNP references the reference folder
						if skip_mauve parameter is True this the align function will only format xmfa file paths
					'''
					xmfa_files = self.align(q)
					LOGGER.debug(xmfa_files)
					if len(xmfa_files) == 0: ## if keep going is set and mauve exits with an error continue to next sequence
						LOGGER.debug("Mauve exited with a non zero exit status, continue with next sample!")
						LOGGER.warning("Mauve error skip {sample}".format(q))
						self.xmfa_files = []
						continue
					'''Parse Mauve XMFA output and find SNPs; returns SNPS (for the visual tree) and SNP_info (text file output)'''
					LOGGER.info("Find SNPs")
					try:
						SNPS,SNP_info,called_snps = self.find_snps_multiproc(xmfa_obj=parse_xmfa_obj,xmfa_files=xmfa_files,export=True)
					except FileNotFoundError:
						LOGGER.warning("One or several xmfa files were not found for {qfile} continue with next file".format(qfile=qfile))
						self.xmfa_files = []
						continue
					'''If file export is requested print the result for each SNP location to file'''
					if self.export:
						outputfile = "{outdir}/{xmfa}_not_called.txt".format(outdir=self.outdir,xmfa=self.query_name)
						outputfile2 = "{outdir}/{xmfa}_snps.txt".format(outdir=self.outdir,xmfa=self.query_name)

						LOGGER.info("Printing SNP info of non called SNPs to {file}".format(file=outputfile))
						self.csnpdict = {}
						'''Print SNPs to tab separated file'''
						with open(outputfile,"w") as snplist_out:
							print("\t".join(["Name","Reference","Pos","Ancestral base","Derived base", "Target base"]),file=snplist_out)
							for snp in SNP_info:
								if snp[0] in called_snps:
									self.csnpdict[snp[0]] = snp
								else:
									print("\t".join(snp),file=snplist_out)

					'''If save tree is requested print tree using ETE3 prints a pdf tree output'''
					SNP = "NA" ## Default message if SNP cannot be confirmed
					final_snp,message,called = self.create_tree(SNPS,self.query_name,called_snps,self.save_tree,min_required_hits=self.min_required_hits,strictness=self.strictness)
					if final_snp:
						SNP = final_snp[1]
						if not final_snp[2][0]:  ## if snp was never confirmed print NA
							SNP = "NA"
						if self.export:
							with open(outputfile2, "w") as called_out:
								if True:
									print("\t".join(["Name","Reference","Pos","Ancestral base","Derived base", "Target base"]),file=called_out)
									for snp in called:
										print("\t".join(self.csnpdict[snp[1]]),file=called_out)
								print("SNP path: {path}".format(path=";".join([snp[1] for snp in called])),file=called_out)
								print("Final SNP: {snp} found/depth: {found}/{depth}".format(snp=SNP,depth=int(final_snp[0]),found=final_snp[2][1]),file=called_out)
						LOGGER.info("Final SNP: {snp} found/depth: {found}/{depth}".format(snp=SNP,depth=int(final_snp[0]),found=final_snp[2][1]))
					else:
						if self.export:
							with open(outputfile2, "a") as called_out:
								print("Final SNP: {snp}".format(snp=SNP), file=called_out)
						LOGGER.info(message)
					if self.summary and SNP != "NA":
						self.summary_set |= set([SNP])
						self.called_genome[SNP] = self.query_name
					if self.export:
						print("{query}: {SNP}".format(query=self.query_name, SNP=SNP))
					'''Clean references to aligned xmfa files between queries if several was supplied'''
					self.xmfa_files = []
				except Exception as e:
					if not self.keep_going:
						LOGGER.warning("An error occured during processing of {file}: {exception}".format(file=self.query_name, exception=e))
						raise CanSNPer2Error("A file did not run correctly exit CanSNPer2 (use --keep_going to continue with next file!)")
					else:
						LOGGER.debug("An error occured during processing of {file}: {exception}".format(file=self.query_name, exception=e))

		if self.summary:
			self.print_summary()
		'''Finally clean up temporary folder when all alignments and trees has been printed!'''
		if not self.keep_temp and len(self.query) > 0: ## if keep temp is turned on do not remove away alignments also if no input files were given
			self.cleanup()

		LOGGER.info("CanSNPer2 finished successfully, files can be found in {outdir}".format(outdir=self.outdir+"/"))
