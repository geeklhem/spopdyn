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

## 2 SP. equal sigma, equal muH. 
param["name"] = "EX03_01"
species.append(np.array([(0.5,0.5,0.1,0.1,0),(0.5,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
dt.append((species[0][1,1]**2-species[0][0,1]**2)/(2.*(species[0][1,1]-species[0][0,1])) - temperature[0,0])
pts.append(([0.5],[0.5]))

## N SP. equal sigma, equal muH. 
param["name"] = "EX03_02a"
param["T_range"] = np.linspace(0,0.5,100)
species.append(np.array([(0.5,muT,0.05,0.05,0) for muT in (.1,.2,.3,.4,.5,.6,.7,.8,.9)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5],[0.5]))
dt.append(0)

## N SP. equal sigma, equal muH. 
param["name"] = "EX03_02"
param["T_range"] = np.linspace(0,0.5,100)
species.append(np.array([(0.5,muT,0.5,0.5,0) for muT in (.1,.2,.3,.4,.5,.6,.7,.8,.9)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5],[0.5]))
dt.append(0)

## 4 SP. equal sigma, diff muH. 
param["name"] = "EX03_03"
param["T_range"] = np.linspace(0,0.4,80)
habitat = np.zeros((2,2))+0.4
habitat[:,0] = 0.6
species.append(np.array([(0.4,0.5,0.1,0.1,0),
                         (0.6,0.4,0.1,0.1,0),
                         (0.4,0.7,0.1,0.1,0),
                         (0.6,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5,0.5],[0.4,0.6]))
dt.append(0)

## 2 SP. equal sigma, diff muH. 
param["name"] = "EX03_04"
param["T_range"] = np.linspace(0,0.4,80)
habitat = np.zeros((2,2))+0.5
habitat[:,0] = 0.6
species.append(np.array([(0.5,0.5,0.1,0.1,0),
                         (0.6,0.8,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5,0.5],[0.5,0.6]))
dt.append(0)

## 2 SP. equal sigma, diff muH. 
param["name"] = "EX03_04b"
param["T_range"] = np.linspace(0,0.4,80)
habitat = np.zeros((2,2))+0.5
habitat[:,0] = 0.6
species.append(np.array([(0.5,0.5,0.1,0.1,0),
                         (0.6,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5,0.5],[0.5,0.6]))
dt.append(0)


## 2 SP. equal sigma, diff muH. 
param["name"] = "EX03_04c"
param["T_range"] = np.linspace(0,0.5,80)
habitat = np.zeros((2,2))+0.5
habitat[:,0] = 0.8
species.append(np.array([(0.5,0.5,0.1,0.1,0),
                         (0.8,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))
pts.append(([0.5,0.5],[0.5,0.8]))
dt.append(0)



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
    plt.vlines(dti,*plt.ylim())
plt.savefig("ex03.png")
plt.savefig("ex03.eps")
