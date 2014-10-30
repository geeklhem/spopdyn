"""
DCTI = F(DT) for only one habitat
"""

from __future__ import division
import numpy as np
import spopdyn.display as d
from spopdyn.dtdcti import applyDT 
import matplotlib.pyplot as plt

temperature = np.zeros((2,2))+0.5
habitat = np.zeros((2,2))+0.5

param = {"alpha":1,
         "d":0,
         "K":500,
         "tend":25,
         "dt":0.1,
         "replicas":1,
         "frames":250,
         "m":1e-1,
         "T_range":np.linspace(0,0.2,40),
}

species = []
outputs = []
dt = []
pts = []

def tc_2sp(T,muT1,muT2):
    return (muT1 + muT2)/2 - T

## 2 SP. equal sigma, equal muH. 
param["name"] = "EX03_01"
species.append(np.array([(0.5,0.5,0.1,0.1,0),(0.5,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
dt.append(tc_2sp(temperature[0,0],species[-1][0,1],species[-1][1,1]))
pts.append(([0.5],[0.5]))

## 2 SP. equal sigma, equal muH. 
param["name"] = "EX03_01b"
species.append(np.array([(0.5,0.5,0.1,0.1,0),(0.5,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
dt.append(tc_2sp(temperature[0,0],species[-1][0,1],species[-1][1,1]))
pts.append(([0.5],[0.5]))

## N SP. equal sigma, equal muH. 
param["name"] = "EX03_02a"
param["T_range"] = np.linspace(0,0.5,100)
possible_mut = (.1,.2,.3,.4,.5,.6,.7,.8,.9)
species.append(np.array([(0.5,muT,0.05,0.05,0) for muT in possible_mut]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5],[0.5]))
dt.append([tc_2sp(temperature[0,0],species[-1][x,1],species[-1][y,1])
           for x,y in zip(range(len(possible_mut))[:-1],range(len(possible_mut))[1:])
              ])

## N SP. equal sigma, equal muH. 
param["name"] = "EX03_02"
param["T_range"] = np.linspace(0,0.5,100)
species.append(np.array([(0.5,muT,0.5,0.5,0) for muT in possible_mut]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5],[0.5]))
dt.append([tc_2sp(temperature[0,0],species[-1][x,1],species[-1][y,1])
           for x,y in zip(range(len(possible_mut))[:-1],range(len(possible_mut))[1:])
              ])



color = "rb"
N = len(outputs)
plt.figure(figsize=(10,N*5))

for n,(out,sp,dti,pt) in enumerate(zip(outputs,species,dt,pts)):
    plt.subplot(N,2,2*n+1)
    d.niches(sp)
    plt.scatter(*pt)
    plt.subplot(N,2,2*(n+1))
    plt.xlabel("$\Delta T$")
    plt.ylabel("$\Delta CTI$")
    plt.scatter(out["deltaT"],out["deltaCTI"],c=[color[x<0] for x in out["deltaT"]])
    plt.hlines(0,*plt.xlim())
    plt.vlines(0,*plt.ylim())
    plt.vlines(dti,*plt.ylim(),color="grey")
plt.savefig("ex03.png")
plt.savefig("ex03.eps")
