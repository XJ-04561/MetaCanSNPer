

def createLogger(name):
    from MetaCanSNPer.Globals import LOGGER
    return LOGGER.getChild(name.split(".")[-1])