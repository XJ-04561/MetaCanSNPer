

from functools import cached_property, cache
import sqlite3, hashlib, re
from typing import Generator, Callable, Iterable, Self, overload, final, Literal, Any

whitespacePattern = re.compile("\s+")

from MetaCanSNPer.Globals import *
from MetaCanSNPer.modules.Databases._Constants import *
from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger("Databases")

DATABASE_VERSION_HASH : str = ""
STRICT : bool = False