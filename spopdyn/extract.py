from __future__ import division 
import logging
import numpy as np
import subprocess
logger = logging.getLogger("spopdyn")

def data_extract_single(fi):
    abundance = []
    
    with open(fi) as f:
        for l in f:
            if len(l) > 1:
                abundance.append(l.split(","))
    abundance = np.array(abundance)
    plt.plot(abundance)
    plt.xlabel("Time")
    plt.ylabel("Abundance")
    plt.tight_layout()
    plt.show()
    plt.clf()

def data_extract_2D(fi,n):
    biomass = []
    nbspecies = []
    diversity = []
    t = -1

    # remove empty lines
    subprocess.check_call("sed -i '/^$/d' {}".format(fi),shell=True)
    
    with open(fi) as f:
        for i, l in enumerate(f):
            pass
    Tmax = i+1
    tt = int(Tmax / n)

    
    with open(fi) as f:
        for t,l in enumerate(f):
            if t % tt == 0:

                loc = l.split("\t")

                gridpoints = np.sqrt(len(loc))

                species = len(loc[0].split(","))

                biomass.append(np.zeros((gridpoints,gridpoints),
                                          dtype=int))
                nbspecies.append(np.zeros((gridpoints,gridpoints),
                                          dtype=int))
                diversity.append(np.zeros((gridpoints,gridpoints),
                                          dtype=float))

                for j,ab_by_sp in enumerate(loc):
                    ab_by_sp = np.fromstring(ab_by_sp,sep=",",dtype=int)
                    ab_by_sp = ab_by_sp[ab_by_sp.nonzero()]

                    biomass[-1].flat[j] = ab_by_sp.sum()
                    nbspecies[-1].flat[j] = len(ab_by_sp)

                    pi = ab_by_sp/ab_by_sp.sum()
                    diversity[-1].flat[j] = -np.sum( pi * np.log(pi))

    return biomass,nbspecies,diversity


