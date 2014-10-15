import cPickle as pkle
import os
import multiprocessing as mp
# Fix pool size
try:
    p = {24:6,8:5}[mp.cpu_count()]
except KeyError:
    p = 1


import matplotlib.pyplot as plt
import numpy as np

import spopdyn.experiment as e
import spopdyn.extract as ext

e.logger.setLevel(e.logging.CRITICAL)

def measure(args):
    dt,CTI,CGI,habitat,temperature,species,initial_pop,param = args
    param["name"] = param["path"] + "/DT{}".format(dt)
    new_temp = temperature+dt
    if not os.path.exists(param["name"]+"/data.pkle"):
        e.experiment(habitat,new_temp,species,initial_pop,param)
        plt.close()
    with open(param["name"]+"/data.pkle","r") as f:
        data = pkle.load(f) 
    return dt,data[0]["Local CTI"][-1]-CTI, data[0]["Local CSI"][-1]-CGI

def applyDT(habitat,temperature,species,param):
    param["path"] = param["name"] + "_applyDT"
    if not os.path.exists( param["path"]):
        os.mkdir( param["path"])
    param["name"] =  param["path"]+"/intial_conditions"
    uniform_pop = [np.zeros(temperature.shape) + 10]*len(species)
    
    e.experiment(habitat,temperature,species,uniform_pop,param)
    
    initial_pop = ext.get_final_pop("{}/PDM/PSSA_trajectory_0_0.txt".format(param["name"]),
                                    param["frames"])
    
    with open("{}/data.pkle".format(param["name"]),"r") as f:
        data = pkle.load(f)
        
        stdcti = np.std(data[0]["Local CTI"][-10:])
        stdcgi = np.std(data[0]["Local CSI"][-10:])
        if stdcgi > 1e-3 or stdcti > 1e-3:
            print "Warning: Initial conditions does not seems to be at the equilibrium"
            print "STD on 10 last timesteps: {} (CTI) and {} (CGI)".format(stdcti,stdcgi)
        CTI  = data[0]["Local CTI"][-1]
        CGI  = data[0]["Local CSI"][-1]
        
    # Measure the DCTI,DCGI for different DT.
    args = []
    if "T_range" not in param:
        param["T_range"] = np.linspace(0.05,0.5,10)
        
    for T in param["T_range"]:
        for sign in (1,-1):
            dt = sign*T
            args.append((dt, CTI, CGI, habitat, temperature, species, initial_pop, param))
    DT,DCTI,DCGI = zip(*POOL.map(measure,args))
    
    out = {"deltaT":DT,
            "deltaCTI":DCTI, "deltaCGI":DCGI,
            "initCTI":CTI, "initCGI":CGI}
    
    with open( param["path"]+"/dt_dcti_dcgi.pkle","w") as f:
        pkle.dump(out,f) 
    return out

POOL = mp.Pool(p)
