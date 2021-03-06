"""
DCTI = F(T) for two different habitats. 
"""
from __future__ import division
import numpy as np
import spopdyn.display as d
from spopdyn.dtdcti import applyDT
from spopdyn.provenance import metadata
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
tc = []

def tc_2sp(muT1,muT2):
    return (muT1 + muT2)/2 


def tc_2habs(muH1,muH2,muT1,muT2):
    """Return Tc for H=muH1 et H=muH2"""
    tc_muH1 = (muH1 - muH2)**2 / (2*(muT2 - muT1))
    tc_muH2 = (muH2 - muH1)**2 / (2*(muT1 - muT2))
    tc_muH1 += tc_2sp(muT1,muT2)
    tc_muH2 += tc_2sp(muT1,muT2)
    return (tc_muH1,tc_muH2)

## 2 SP. equal sigma, diff muH. 
param["name"] = "EX04_01"
param["T_range"] = np.linspace(0,0.4,80)
habitat = np.zeros((2,2))+0.5
habitat[:,0] = 0.6
species.append(np.array([(0.5,0.5,0.1,0.1,0),
                         (0.6,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5,0.5],[0.5,0.6]))
tc.append(tc_2habs(species[-1][0,0], species[-1][1,0],
                   species[-1][0,1], species[-1][1,1]))
dt.append([x-temperature[0,0] for x in tc[-1]])

## 2 SP. equal sigma, diff muH. 
param["name"] = "EX04_02"
param["T_range"] = np.linspace(0,0.5,80)
habitat = np.zeros((2,2))+0.5
habitat[:,0] = 0.6
species.append(np.array([(0.5,0.5,0.1,0.1,0),
                         (0.6,0.8,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5,0.5],[0.5,0.6]))
tc.append(tc_2habs(species[-1][0,0], species[-1][1,0],
                   species[-1][0,1], species[-1][1,1]))
dt.append([x-temperature[0,0] for x in tc[-1]])

## 4 SP. equal sigma, diff muH. 
param["name"] = "EX04_03"
param["T_range"] = np.linspace(0,0.4,80)
habitat = np.zeros((2,2))+0.4
habitat[:,0] = 0.6
species.append(np.array([(0.4,0.5,0.1,0.1,0),
                         (0.6,0.4,0.1,0.1,0),
                         (0.4,0.7,0.1,0.1,0),
                         (0.6,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5,0.5],[0.4,0.6]))
tc.append((tc_2sp(species[-1][0,1],species[-1][2,1]),
           tc_2sp(species[-1][1,1],species[-1][3,1])))
dt.append([x-temperature[0,0] for x in tc[-1]])

color = "rb"
N = len(outputs)
plt.figure(figsize=(10,N*5))

for n,(out,sp,tci,dti,pt) in enumerate(zip(outputs,species,tc,dt,pts)):
    plt.subplot(N,2,2*n+1)
    d.niches(sp)
    plt.scatter(*pt)
    plt.vlines(tci,*plt.ylim(),linestyle=":",color="darkblue")
    plt.subplot(N,2,2*(n+1))
    plt.xlabel("$\Delta T$")
    plt.ylabel("$\Delta CTI$")
    plt.scatter(out["deltaT"],out["deltaCTI"],c=[color[x<0] for x in out["deltaT"]])
    plt.hlines(0,*plt.xlim())
    plt.vlines(0,*plt.ylim())
    plt.vlines(dti,*plt.ylim(),linestyle=":",color="darkblue")
plt.tight_layout()

plt.savefig("dt_dcti_2habs.png")
plt.savefig("dt_dcti_2habs.eps")
metadata("dt_dcti_2habs.png")
metadata("dt_dcti_2habs.eps")
