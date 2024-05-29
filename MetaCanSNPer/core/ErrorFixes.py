
from MetaCanSNPer.Globals import *
from MetaCanSNPer.core.Wrappers import *

from copy import copy as copyObject
from shutil import copy

class SolutionContainer(Logged):

    def __init__(self, obj : "ProcessWrapper"):
         self.obj = obj

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

    def _11(self):
        outDir = self.obj.Lib.tmpDir.create("progressiveMauve").create("_11")

        self.LOG.info(f"Fixing exitcode 11 by replacing occurrances of '-' with 'N' from {self.obj.Lib.query!r} into separate new files.")
        out = []
        for q in self.obj.Lib.query:
            tmpName = "{}{}".format(os.path.splitext(os.path.basename(q))[0], os.path.splitext(os.path.basename(q))[1])
            self.LOG.debug(f"copy('{q}', '{tmpName}')")
            copy(q, tmpName)
            self.LOG.debug(f"open('{tmpName}', 'r+b')")
            with open(tmpName, "r+b") as f:
                c = b" "
                while c != b"":
                    c = f.read(1)
                    if c == b"-":
                        f.seek(-1, 1)
                        f.write(b"N")
                f.close()
            out.append(tmpName)
        
        self.LOG.info(f"New query: {out}")
        self.obj.Lib = copyObject(self.obj.Lib)
        object.__setattr__(self.obj.Lib, "_lib", object.__getattribute__(self.obj.Lib, "_lib").copy())
        self.obj.Lib.setQuery(out, abs=True)

class minimap2(SolutionContainer):

    def _11(self):
        self.LOG.info("Fixing exitcode 11 by changing paths to symbolic links into real paths.")
        from PseudoPathy import Path
        for i in range(len(self.obj.Lib.query)):
            msg = f"{self.obj.Lib.query[i]!r} now becomes "
            self.obj.Lib.query[i] = Path(os.path.realpath(self.obj.Lib.query[i]))
            msg += f"{self.obj.Lib.query[i]!r}"
            self.LOG.debug(msg)
        
        for genome, path in self.obj.Lib.references.items():
            msg = f"{path!r} now becomes "
            self.obj.Lib.references[genome] = Path(os.path.realpath(path))
            msg += f"{path!r}"
            self.LOG.debug(msg)

def get(softwareName) -> SolutionContainer:
	for c in SolutionContainer.__subclasses__():
		if c.__name__ == softwareName:
			return c
	return SolutionContainer

try:
    from MetaCanSNPer.core.Wrappers import ProcessWrapper
except ImportError:
    pass