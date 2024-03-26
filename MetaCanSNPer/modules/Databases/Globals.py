

from functools import cached_property
import sqlite3
from typing import Generator, Callable, Iterable, Self, overload, final, Literal, Any

from MetaCanSNPer.Globals import *
from MetaCanSNPer.modules.Databases._Constants import *
from MetaCanSNPer.modules.LogKeeper import createLogger

LOGGER = createLogger("Databases")