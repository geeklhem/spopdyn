import logging
import subprocess
import os 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


logger = logging.getLogger("spopdyn")

def exp_summary(habitat,temperature,species):
    plt.subplot(2,2,1)
    niches(species)
    plt.subplot(2,2,2)
    environment(habitat,temperature)
    plt.subplot(2,2,3)
    plt.imshow(habitat,
               interpolation='None',
               cmap="Greens",
               vmin=0,vmax=1,
               aspect="equal")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Habitat quality")
    plt.colorbar(fraction=0.05)
    plt.subplot(2,2,4)
    plt.imshow(temperature,
               interpolation='None',
               cmap="bwr",
               vmin=0,
               vmax=1,
               aspect="equal",)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.colorbar(fraction=0.05)
    plt.title("Temperature")


def environment(habitat,temperature):
    gp = habitat.shape[0]
    img = np.zeros((gp,gp,3))
    img[:,:,2] = np.floor(255*(temperature))
    img[:,:,0] = np.floor(255*(1-temperature))
    img[:,:,1] = np.floor(125*(1-habitat))
    plt.imshow(img,
               aspect="equal",
               interpolation='None')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Environmental value")

def niches(species,gp=25):
    img = np.zeros((gp,gp,3))

    for x in range(gp):
        img[-x-1,:,1] = (x+1)*int(125/gp)
        img[:,-x-1,0] = (x+1)*int(255/gp)
        img[:,x,2] = (x+1)*int(255/gp)

    plt.imshow(img,
               interpolation='None',
               extent=(0,1,1,0),
               aspect="equal",)
    ax = plt.gca()
    for sp in species:
        ellipse = Ellipse(xy=(sp[0], sp[1]), width=sp[2], height=sp[3], 
                          edgecolor='k',alpha=0.1, fc='None', lw=1,angle=sp[4]*360)
        ax.add_patch(ellipse)

    plt.xlabel("Temperature")
    plt.ylabel("Habitat quality")
    plt.title("Ecological niches")


def matrix_movie(matrices,titles,maxvalues,path):
    cmap_list = ["Blues","Greens","Oranges","RdPu","Reds"]
    colors = ["blue","green","orange","pink","red"]
    logger.info("Images...")
    Tmax = len(matrices[0])
    n = len(matrices)
    means = np.zeros((n,Tmax))
    var =  np.zeros((n,Tmax))
    
    
    for t in range(Tmax):
        print "{}/{} {}[A".format(t,Tmax,chr(27))
        if not os.path.exists("{}{:05}.png".format(path,t)):
            fig = plt.figure(figsize=(15,5))
            for m,(matrix,title,mx) in enumerate(zip(matrices,titles,maxvalues)):
                plt.subplot(2,n,m+1)
                plt.imshow(matrix[t], aspect="equal",
                           vmin=0, vmax=mx,
                           cmap=cmap_list[m],
                           interpolation="None")
                plt.colorbar(fraction=0.05)
                plt.title("{}".format(title))

                plt.subplot(3,n,n+m+1)
                plt.title("{}".format(title))
                means[m,t] = matrix[t].mean()
                plt.plot(means[m,:t],color=colors[m])
                plt.xlim(0,Tmax)
                
                plt.subplot(3,n,n+2*m+1)
                var[m,t] = matrix[t].var()
                plt.plot(var[m,:t],"-.",color=colors[m])
                plt.xlim(0,Tmax)

            plt.tight_layout()
            plt.savefig("{}{:05}.png".format(path,t),pad_inches=0)
            plt.clf() 

    logger.info("Movie...")
    if not os.path.exists("{0}video.mp4".format(path)):
        subprocess.call("avconv -y -r 3 -i {0}%05d.png {0}video.mp4".format(path), shell=True)
    subprocess.call("mplayer {}video.mp4  -idle -fixed-vo".format(path), shell=True)


