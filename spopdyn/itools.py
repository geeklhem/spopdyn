""" itools.py - a suite of interactive tools for rapid analysis of spopdyn output """
from __future__ import division
import logging
import subprocess
import os 

import numpy as np
import matplotlib.pyplot as plt
import spopdyn.extract as extract

logger = logging.getLogger("spopdyn")

def popdyn(fi):
    """ Extract population dynamics from a libpssa file and display it """
    p = extract.popdyn(fi)
    plt.xlabel("Time")
    plt.ylabel("Abundance")
    plt.plot(p)
    plt.show()
