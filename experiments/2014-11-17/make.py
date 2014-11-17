"""
Sensitivity analysis. Usibng Morris Method
"""
from __future__ import division
import os
import cPickle as pkle
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt

from SALib.sample import morris_oat
from SALib.analyze import morris

from spopdyn.provenance import metadata
from spopdyn.landscape import seq_gaussian
import spopdyn.experiment as e

def measure_tau(T,rho,R,D,s):

    param = {
        "alpha":1,
        "K":500,
        "tend":25,
        "dt":0.1,
        "replicas":1,
        "frames":250,
        "m":1e-1,
        "name":"SA_",
        "step_T":.25
    }
    
    #---- Folder ----#
    param["path"] = param["name"] + hash((T,rho,R,D,s))
    if not os.path.exists( param["path"]):
        os.mkdir( param["path"])

    #---- Parameters ----#
    temperature = np.zeros((10,10))+T
    habitat = seq_gaussian(10,rho)
    species = np.zeros((R,5))
    species[:,:2] = np.random.random(size=(R,2))
    species[:,2:4] = np.random.random(size=(R,2))*s 
    param["d"] = D

    #---- Go to equilibrium -----#
    param["name"] =  param["path"]+"/intial_conditions"
    _,initial_pop = e.equilibrium(habitat,temperature,species,None,param)
    
    #---- Do the thermal step ----#
    param["name"] =  param["path"]+"/step"
    new_temp = temperature+param["step_T"]
    data,_ = e.equilibrium(habitat,new_temp,species,initial_pop,param)
    CTI  = data[0]["Local CTI"]
    for i in range(len(CTI))[:-20]:
        if np.std(CTI[i:i+20])<1e-3:
            tau = i 
            break
    #---- Measure tau ----# 
    return tau 


def main():
    #print param
    param = np.array([
        ["T",0,.75],
        ["rho",0,1],
        ["R",10,50],
        ["D",0,5],
        ["s",0,.6],
    ])
    np.savetxt('parameters.txt',param,fmt='%s',delimiter=' ')
    # Generate samples
    param_values = morris_oat.sample(100, "parameters.txt", num_levels=10, grid_jump=5)
    # Save the parameter values in a file (they are needed in the analysis)
    np.savetxt('model_input.txt', param_values, delimiter=' ')

    
    Y = np.array(POOL.map(measure_tau,param_values))
    np.savetxt("model_output.txt", Y, delimiter=' ')

    # Perform the sensitivity analysis using the model output
    Si = morris.analyze('parameters.txt', 'model_input.txt', 'model_output.txt',
                        column=0, conf_level=0.95, print_to_console=False)

    with open("data.pkle","w") as f:
        pkle.dump({"X": param_values, "Y":Y, "Morris":Si},f)
    metadata("data.pkle")

    
    plt.scatter(Si["mu_star"],Si["sigma"])

    for x,y,s in zip(Si["mu_star"],Si["sigma"],param[:,0]):
        
        plt.annotate(s,(x,y), xytext=(5, 5),textcoords='offset points')
            
    plt.xlabel(r"$\mu^*$")
    plt.ylabel(r"$\sigma$")
        
    plt.savefig("sensitivity_with_morris.png")
    metadata("sensitivity_with_morris.png")

    plt.savefig("sensitivity_with_morris.eps")
    metadata("sensitivity_with_morris.eps")
        
    return Si


# Fix pool size
try:
    p = {24:6,8:5}[mp.cpu_count()]
except KeyError:
    p = 1
POOL = mp.Pool(p)
if __name__ == "__main__":
    main()

