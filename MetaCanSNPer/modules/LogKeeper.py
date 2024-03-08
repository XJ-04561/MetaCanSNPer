
from MetaCanSNPer.Globals import *


def createLogger(name):
    logger = logging.Logger("MetaCanSNPer."+name)
    logger.addHandler(LOGGER_FILEHANDLER)

    return logger