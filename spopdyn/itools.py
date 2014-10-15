from __future__ import division
import logging
import subprocess
import os 

import numpy as np
import matplotlib.pyplot as plt
import spopdyn.extract as extract

logger = logging.getLogger("spopdyn")

def popdyn(fi):
    p = extract.popdyn(fi)
    plt.plot(p)
    plt.show()
