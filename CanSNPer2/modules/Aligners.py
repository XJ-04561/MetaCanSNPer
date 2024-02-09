

import logging
import LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)
from Wrappers import Aligner

'''
	All that is needed to create a new implementation is to inherit from the correct software type ('Aligner' in this case) and set
	the two class attributes accordingly. See 'Timeout' and 'Sleep' for a minimalist example.
'''

class Timeout(Aligner):
	softwareName = "timeout"
	commandTemplate = "timeout {ref} {target}"

class Sleep(Aligner):
	softwareName = "sleep"
	commandTemplate = "sleep {ref} {target}"

class ProgressiveMauve(Aligner):
	softwareName = "progressiveMauve"
	commandTemplate = "progressiveMauve {ref} {target}"

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