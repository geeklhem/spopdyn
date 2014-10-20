import numpy as np
import spopdyn.display as d
from spopdyn.dtdcti import applyDT 
import spopdyn.sensitivity as sty
import matplotlib.pyplot as plt


temperature = np.zeros((2,2))+0.5
habitat = np.zeros((2,2))+0.5

param = {"alpha":1,
         "name":"ex02",
          "d":1,
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
species.append(np.array([(0.5,0.5,0.1,0.1,0),(0.5,0.6,0.1,0.1,0)]))
outputs.append(applyDT(habitat,temperature,species[-1],param))

outputs[-1]["S"] = sty.sensitivity_analysis(np.array(outputs[-1]["deltaT"]),
                                            np.array(outputs[-1]["deltaCTI"]))
dt = (species[0][1,1]**2-species[0][0,1]**2)/(2.*(species[0][1,1]-species[0][0,1])) - temperature[0,0]

color = "rb"
N = len(outputs)
plt.figure(figsize=(10,N*5))

for n,(out,sp) in enumerate(zip(outputs,species)):
    plt.subplot(N,2,2*n+1)
    d.niches(sp)
    plt.subplot(N,2,2*(n+1))
    sty.plot_poly(out["S"])
    plt.xlabel("$\Delta T$")
    plt.ylabel("$\Delta CTI$")
    plt.scatter(out["deltaT"],out["deltaCTI"],c=[color[x<0] for x in out["deltaT"]])
    plt.hlines(0,*plt.xlim())
    plt.vlines(0,*plt.ylim())
    plt.vlines(dt,*plt.ylim())
plt.savefig("ex_02_m.png")
plt.savefig("ex_02_m.eps")
