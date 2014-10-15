import numpy as np
import spopdyn.display as d
from spopdyn.dtdcti import applyDT 
import matplotlib.pyplot as plt

temperature = np.zeros((2,2))+0.5
habitat = np.zeros((2,2))+0.5
species = np.array([(0.5,0.5,0.1,0.1,0),(0.6,0.5,0.1,0.1,0)])
param = {"alpha":1,
         "name":"2sp",
          "d":1,
          "K":500,
          "tend":25,
          "dt":0.1,
          "replicas":1,
          "frames":250,
          "m":1e-1,
          "T_range":np.linspace(0,0.2,10),
          }
dt = (species[1,0]**2-species[0,0]**2)/(2.*(species[1,0]-species[0,0])) - temperature[0,0]
out = applyDT(habitat,temperature,species,param)
color = "rb"

plt.figure(figsize=(10,5))
plt.subplot(121)
d.niches(species)
plt.subplot(122)
plt.xlabel("$\Delta T$")
plt.ylabel("$\Delta CTI$")
plt.scatter(out["deltaT"],out["deltaCTI"],c=[color[x<0] for x in out["deltaT"]])
plt.hlines(0,*plt.xlim())
plt.vlines(0,*plt.ylim())
plt.vlines(dt,*plt.ylim())
plt.savefig("ex01.png")
plt.savefig("ex01.eps")
