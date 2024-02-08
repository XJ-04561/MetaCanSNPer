from Wrappers import *
import os.path

import logging
import CanSNPer2.modules.LogKeeper as LogKeeper

LOGGER = LogKeeper.createLogger(__name__)

class SolutionContainer:
    def __contains__(self, key):
        try:
            self[key]
            return True
        except:
            return False

    def __getitem__(self, key):
        if type(key) is not int :
            raise TypeError("Container for solutions can only be indexed with type 'int'")
        
        try:
            ret = self.__getattribute__("_{sign}{number}".format(number=abs(key), sign="_" if key<0 else ""))
        except:
            raise KeyError("No solution available for returncode {}.".format(key))
        return ret


class progressiveMauve(SolutionContainer):
    @staticmethod
    def _11(obj : ProgressiveMauve, offenders : list[int]):
        LOGGER.debug("sed 's/-/N/g' {query} > {tmpName}.tmp".format(query=obj.query, query_base=obj.queryName, tmpdir=obj.tmpdir, sep=os.path.sep))
        tmpName = "{}.tmp".format(os.path.join([obj.Lib.tmpDir, obj.queryName]))
        os.system("sed 's/-/N/g' {query} > {tmpName}".format(query=obj.Lib.getQuery(), tmpName=tmpName))
        
        LOGGER.debug("New query: {nq}".format(nq=tmpName))
        obj.Lib.setQuery(tmpName, abs=True)
    

Indexers = {
    "progressiveMauve" : progressiveMauve
}