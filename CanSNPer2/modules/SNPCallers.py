

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
	commandTemplate = "timeout {ref} {target}"

class Sleep(SNPCaller):
	softwareName = "sleep"
	commandTemplate = "sleep {ref} {target}"


def get(softwareName) -> SNPCaller:
	for c in SNPCaller.__subclasses__():
		if c.softwareName == softwareName:
			return c
	return None