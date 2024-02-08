import logging

def createLogger(name):
    fileHandler = logging.FileHandler("MetaCanSNPer.log")
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.setFormatter(logging.Formatter("[%(name)s] %(asctime)s - %(levelname)s: %(message)s"))
    logger = logging.Logger("MetaCanSNPer."+name)
    logger.addHandler(fileHandler)

    return logger