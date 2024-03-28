

from typeguard import install_import_hook

install_import_hook("MetaCanSNPer.modules.Databases.Globals")
import MetaCanSNPer.modules.Databases.Globals as Globals

install_import_hook("MetaCanSNPer.modules.Databases.Databases")
from MetaCanSNPer.modules.Databases.Databases import openDatabase, DatabaseReader, DatabaseWriter

install_import_hook("MetaCanSNPer.modules.Databases.Columns")
from MetaCanSNPer.modules.Databases.Columns import ColumnFlag
import MetaCanSNPer.modules.Databases.Columns as Columns

install_import_hook("MetaCanSNPer.modules.Databases.Functions")
from MetaCanSNPer.modules.Databases.Functions import downloadDatabase

install_import_hook("MetaCanSNPer.modules.Databases.Tables")
from MetaCanSNPer.modules.Databases.Tables import SNPTable, ReferenceTable, NodeTable, TreeTable

install_import_hook("MetaCanSNPer.modules.Databases.Tree")
from MetaCanSNPer.modules.Databases.Tree import Branch

install_import_hook("MetaCanSNPer.modules.Databases.Test")
import MetaCanSNPer.modules.Databases.Test as Test
