from __future__ import division
import logging
import subprocess
import os 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA


logger = logging.getLogger("spopdyn")

def show_matrix(matrix,kind="temperature"):
    """
    Display the matrix using matplotlib, the title and the colormap of
    the plot are deduced from the 'kind' parameter.

    Args:
        matrix (np.array): 
        kind (str): kind of information.
    """
    if kind=="temperature":
        cmap = "bwr"
        plt.title("Temperature")
    elif kind=="habitat":
        cmap = "Greens"
        plt.title("Habitat")
    else:
        cmap = "Blues"
    plt.imshow(matrix,
               interpolation='None',
               cmap=cmap,
               vmin=0,
               vmax=1,
               aspect="equal",)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.xticks([])
    plt.yticks([])
    plt.colorbar(orientation="horizontal", fraction=0.045)



def exp_summary(habitat,temperature,species):
    """
    Plot the environmental matrices and a ecological niche representation.

    Args:
        habitat (np.array): habitat value.
        temperature (np.array): temerature value.
        species (np.array): species characteristics for the niche plotting.
    """
    plt.subplot(2,2,1)
    niches(species)
    plt.subplot(2,2,2)
    environment(habitat,temperature)
    plt.subplot(2,2,3)
    show_matrix(habitat,"habitat")
    plt.subplot(2,2,4)
    show_matrix(temperature,"temperature")

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

def niches(species,gp=25,path=None):
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

    names =  len(species) < 6

    for i,sp in enumerate(species):
        ellipse = Ellipse(xy=(sp[1], sp[0]), width=sp[3], height=sp[2], 
                          edgecolor='k',alpha=0.1, fc='None', lw=1,angle=sp[4]*360)
        ax.add_patch(ellipse)
        if names:
            plt.annotate("{}".format(i),(sp[0],sp[1]),size=10)
    plt.scatter(species[:,0],species[:,1],color="k",s=20,marker="x")
    plt.xlim((0,1))
    plt.ylim((0,1))

    if path is not None:
        for c1,c2 in zip(path[:-1],path[1:]):
            plt.arrow(c1[0],c1[1],c2[0]-c1[0],c2[1]-c1[1],color="blue")
    
    plt.xlabel("Temperature")
    plt.ylabel("Habitat")
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


def experimental_report(environment, species, time_series,path=None,events=None):
    """
    Save a report of the experiment as an eps file.

    Input:
        environment (list of tuple): list of environmental variable (name, value matrix),

    """


    M = len(environment)+1
    L = int(np.ceil(1 + len(time_series)/2))
    fig = plt.figure(figsize=(5*M,5*L))
    
    colormaps = ["Greens","bwr","Blues","Oranges","RdPu","Reds"]
    for i,(k,v) in enumerate(environment):
        plt.subplot(L,M,i+1)
        plt.imshow(v,
                   interpolation='None',
                   cmap=colormaps[i%len(colormaps)],
                   vmin=0,vmax=1,
                   aspect="equal")
        plt.xticks([])
        plt.yticks([])
        plt.title(k)
        plt.colorbar(orientation="horizontal", fraction=0.045)
    plt.subplot(L,M,M)
    niches(species,path=path)

    colors = ["blue","green","brown","purple","red"]
    host = [host_subplot(L*100+10+2+j, axes_class=AA.Axes) for j in range(L-1)]


    for i,(k,v) in enumerate(time_series):
        #if False and i%2 != 0:
        #    ax = host[int(i/2)].twinx()
        #else:
        ax = host[int(i/2)]
        ax.set_ylabel(k)
        if len(v) == 2:
            T = len(v[0])
            ax.plot(v[0],
                    label=k,
                    color=colors[i%len(colors)],
                    linewidth=2)
            ax.fill_between(range(len(v[0])),
                              v[0]-v[1], v[0]+v[1],
                              alpha=0.3,
                              color=colors[i%len(colors)])
        else:
            T = len(v)
            ax.plot(range(len(v)),v, color=colors[i%len(colors)], label=k)
            
    
    for h in host:
        h.set_xlim((0,T-1))
        h.legend()
        h.set_xlabel("Time")
        
        h.set_ymargin(0.05)
        h.autoscale(enable=True, axis=u'both', tight=False)

        if events is not None:
            h.vlines(events,*h.get_ylim(),alpha=0.1)


if __name__ == "__main__":
    t_s = [
           ("Richness",  np.random.uniform(size=10)),
           ("Diversity", (np.random.uniform(size=10) , np.random.uniform(size=10)*0.1) ),
           ("CTI",  np.random.uniform(size=10)),
           ("CSI",  np.random.uniform(size=10)),
        ("Biomass", np.random.uniform(size=10)),]
    
    e = [("Habitat", np.random.uniform(0,1,size=(25,25))),
         ("Temperature", np.random.uniform(0,1,size=(25,25)))]

    s = np.random.uniform(0,1,size=(1,5))
    s[:,-1] = 0

    experimental_report(e,s,t_s,".")









