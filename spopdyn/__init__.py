import logging
from spopdyn.itools import popdyn

logger = logging.getLogger("spopdyn")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter('%(levelname)s:%(name)s: %(message)s'))
logger.addHandler(ch)

logging.addLevelName( logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
logging.addLevelName( logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))
logging.addLevelName( logging.INFO, "\033[1;42m%s\033[1;0m" % logging.getLevelName(logging.INFO))
logging.addLevelName( logging.DEBUG, "\033[1;43m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
logging.getLogger().addFilter(logging.Filter("vcontact"))


