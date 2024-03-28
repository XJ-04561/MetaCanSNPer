
from MetaCanSNPer.modules.Databases.Globals import *
import MetaCanSNPer.modules.Databases.Globals as Globals
import MetaCanSNPer.modules.Databases.Columns as Columns
from MetaCanSNPer.modules.Databases.Columns import ColumnFlag
from MetaCanSNPer.modules.Databases._Constants import *

class Branch:

	_connection : sqlite3.Connection
	nodeID : int
	parameters = [
		TREE_COLUMN_CHILD,
		TABLE_NAME_TREE,
		TREE_COLUMN_PARENT
	]

	def __init__(self, connection : sqlite3.Connection, nodeID : int):
		self._connection = connection
		self.nodeID = nodeID
	
	@property
	def children(self) -> Generator[Self,None,None]:
		for (childID,) in self._connection.execute(f"SELECT {TREE_COLUMN_CHILD} FROM {TABLE_NAME_TREE} WHERE {TREE_COLUMN_PARENT} = ?", [self.nodeID]):
			yield Branch(self._connection, childID)