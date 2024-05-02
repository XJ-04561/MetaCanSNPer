

def createLogger(name):
    import logging
    from MetaCanSNPer.Globals import LOGGER_FILEHANDLER
    logger = logging.Logger("MetaCanSNPer."+name)
    logger.addHandler(LOGGER_FILEHANDLER)

    return logger