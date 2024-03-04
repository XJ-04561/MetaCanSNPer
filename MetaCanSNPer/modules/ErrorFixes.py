import os.path, logging
from shutil import copy

class Aligner: pass
class Mapper: pass
class SNPCaller: pass

try:
    import MetaCanSNPer.modules.LogKeeper as LogKeeper
    from MetaCanSNPer.modules.Wrappers import *
except:
    import LogKeeper as LogKeeper
    from Wrappers import *

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
    def _11(obj : Aligner):
        tmpName = "{}.tmp".format(os.path.join([obj.Lib.tmpDir, obj.queryName]))
        LOGGER.info("Fixing exitcode 11 by replacing occurrances of '-' with 'N' from '{query}' into '{tmpName}'.".format(query=obj.Lib.query, tmpName=tmpName))

        LOGGER.debug("copy('{}', '{}')".format(obj.Lib.query, tmpName))
        copy(obj.Lib.query, tmpName)
        LOGGER.debug("open('{}', 'r+b')".format(tmpName))
        with open(tmpName, "r+b") as f:
            c = b" "
            while c != b"":
                c = f.read(1)
                if c == b"-":
                    f.seek(-1, 1)
                    f.write(b"N")
            f.close()
        
        LOGGER.info("New query: {nq}".format(nq=tmpName))
        obj.Lib.setQuery(tmpName, abs=True)

def get(softwareName) -> SolutionContainer:
	for c in SolutionContainer.__subclasses__():
		if c.__name__ == softwareName:
			return c
	return SolutionContainer