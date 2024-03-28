

from functools import cached_property, cache
import sqlite3, hashlib, re
from typing import Generator, Callable, Iterable, Self, overload, final, Literal, Any

from MetaCanSNPer.Globals import *
from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger("Databases")

DATABASE_VERSIONS : dict[str,int] = {

}
STRICT : bool = False