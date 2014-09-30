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

def data_extract_2D(fi,n,param):
    biomass = []
    nbspecies = []
    diversity = []
    CSI = []
    CTI = []
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

                CSI.append(np.zeros((gridpoints,gridpoints),
                                          dtype=float))
                CTI.append(np.zeros((gridpoints,gridpoints),
                                          dtype=float))


                for j,ab_by_sp in enumerate(loc):
                    ab_by_sp = np.fromstring(ab_by_sp,sep=",",dtype=int)
                    CSI[-1].flat[j] = np.sum(ab_by_sp * param[:,1])
                    CTI[-1].flat[j] = np.sum(ab_by_sp * param[:,2])
                    ab_by_sp = ab_by_sp[ab_by_sp.nonzero()]

                    biomass[-1].flat[j] = ab_by_sp.sum()
                    nbspecies[-1].flat[j] = len(ab_by_sp)

                    pi = ab_by_sp/ab_by_sp.sum()
                    diversity[-1].flat[j] = -np.sum( pi * np.log(pi))
                CSI[-1] /= biomass[-1]
                CTI[-1] /= biomass[-1]
    return biomass,nbspecies,diversity,CSI,CTI


def time_serie(martices):
    """
    Arg:
        matrices (list): list of matrices.
    Return:
        (tuple of np.array): mean and spatial variance time series   
    """
    T = len(matrices)
    mean = np.zeros(T)
    var = np.zeros(T)
    for i,M in enumerate(matrices):
        mean[i] = np.mean(M)
        var[i] = np.var(M)
    return (mean,var)

