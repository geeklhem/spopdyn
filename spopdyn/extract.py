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
    return abundance

def get_final_pop(fi,n):
    subprocess.check_call("sed -i '/^$/d' {}".format(fi),shell=True)
    with open(fi) as f:
        for i, l in enumerate(f):
            pass
    Tmax = i+1
    tt = int(Tmax / n)
    ttmax = tt*n

    final_pop = []
    with open(fi) as f:
        for t, l in enumerate(f):
            if t == ttmax:
                loc = l.split("\t")
                gridpoints = np.sqrt(len(loc))
                n_species = len(np.fromstring(loc[0],sep=",",dtype=int))
                for k in range(n_species):
                    final_pop.append(np.zeros((gridpoints,gridpoints)))
                for j,ab_by_sp in enumerate(loc):
                    ab_by_sp = np.fromstring(ab_by_sp,sep=",",dtype=int)
                    for k,ab in enumerate(ab_by_sp):
                        final_pop[k].flat[j] = ab
                    
                
    return final_pop

def data_extract_2D(fi,n,species):
    biomass = []
    nbspecies = []
    diversity = []
    CSI = []
    CTI = []
    total_species = []
    total_div = []

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

                tt_species = np.zeros(species.shape[0])

                for j,ab_by_sp in enumerate(loc):

                    ab_by_sp = np.fromstring(ab_by_sp,sep=",",dtype=int)

                    CSI[-1].flat[j] = np.sum(ab_by_sp * species[:,2])
                    CTI[-1].flat[j] = np.sum(ab_by_sp * species[:,1])
                    ab_by_sp_withoutZero = ab_by_sp[ab_by_sp.nonzero()]

                    biomass[-1].flat[j] = ab_by_sp_withoutZero.sum()
                    nbspecies[-1].flat[j] = len(ab_by_sp_withoutZero)

                    pi = ab_by_sp_withoutZero/ab_by_sp_withoutZero.sum()
                    diversity[-1].flat[j] = -np.sum( pi * np.log(pi))

                    try:
                        tt_species += ab_by_sp
                    except Exception:
                        print "T={},Location: {}\nTotal species: {}\n Species {}".format(j,len(CSI),tt_species,ab_by_sp)
                        raise Exception
                pi =  tt_species[tt_species.nonzero()]/biomass[-1].sum()
                total_div.append(np.array(-np.sum( pi * np.log(pi))))
                total_species.append(np.array(np.sum(tt_species>0)))
                CSI[-1] /= biomass[-1]
                CTI[-1] /= biomass[-1]

    return {"Local Biomass":biomass,
            "Local Species Richness":nbspecies,
            "Local Diversity":diversity,
            "Local CSI":CSI,
            "Local CTI":CTI,
            "Regional Species Richness":total_species,
            "Regional Diversity":total_div,
    }

def extract_multiple_trajectories_to_timeseries(fi,n,param):
    ts_mean = {}
    ts_svar = {}
    #print "{} trajectories files".format(len(fi))
    for f in fi:
        d = data_extract_2D(f,n,param)
        for k in d.keys():
            m,v = time_serie(d[k])
            if k in ts_mean:
                ts_mean[k].append(m)
                ts_svar[k].append(v)
            else:
                ts_mean[k] = [m]
                ts_svar[k] = [v]
            #print "tsmean[{}]: {}\n m: {}".format(k,len(ts_mean[k]),m.shape)
    return_mean = {}
    return_var = {}
    for k in ts_mean.keys():
        arr = np.array(ts_mean[k])
        #print "arr: {}".format(arr.shape)
        return_mean[k] = np.mean(arr,0)
        return_var[k] = np.var(arr,0)
        arr = np.array(ts_svar[k])
        return_mean["spatial variance of"+k] = np.mean(arr,1)
        return_var["spatial variance of"+k] = np.var(arr,1)
    #print "\n mean {}: \n m: {}".format(k,return_mean[k])
    return return_mean,return_var
        
def time_serie(matrices):
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

