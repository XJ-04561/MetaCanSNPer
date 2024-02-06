from Wrappers import *
import os.path

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
        obj.logger.debug("sed 's/-/N/g' {query} > {tmpdir}{sep}{query_base}.tmp".format(query=obj.query, query_base=obj.queryName, tmpdir=obj.tmpdir, sep=os.path.sep))
        os.system("sed 's/-/N/g' {query} > {tmpdir}{sep}{query_base}.tmp".format(query=obj.query, query_base=obj.queryName, tmpdir=obj.tmpdir, sep=os.path.sep))
        
        new_query = "{tmpdir}{sep}{query}.tmp".format(query=obj.query_base, tmpdir=obj.tmpdir, sep=os.path.sep)
        obj.logger.debug("New query: {nq}".format(nq=new_query))
    

Indexers = {
    "progressiveMauve" : progressiveMauve
}