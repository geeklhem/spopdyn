import logging
import subprocess
import os 

import numpy as np
import matplotlib.pyplot as plt


logger = logging.getLogger("spopdyn")


def matrix_movie(matrices,titles,maxvalues,path):
    cmap_list = ["Blues","Greens","Oranges"]
    colors = ["blue","green","orange"]
    logger.info("Images...")
    Tmax = len(matrices[0])
    n = len(matrices)
    means = np.zeros((n,Tmax))
    var =  np.zeros((n,Tmax))

    for t in range(Tmax):
        print "{}/{} {}[A".format(t,Tmax,chr(27))
        if not os.path.exists("{}{:05}.png".format(path,t)):
            for m,(matrix,title,mx) in enumerate(zip(matrices,titles,maxvalues)):
                plt.subplot(2,n,m+1)
                plt.imshow(matrix[t], aspect="equal",
                           vmin=0, vmax=mx,
                           cmap=cmap_list[m],
                           interpolation="None")
                plt.colorbar(fraction=0.05)
                plt.title("{}".format(title))

                plt.subplot(2,n,n+m+1)
                means[m,t] = matrix[t].mean()
                var[m,t] = matrix[t].var()
                plt.plot(means[m,:t],color=colors[m])
                plt.plot(var[m,:t],"-.",color=colors[m])
                plt.title("{}".format(title))
                plt.xlim(0,Tmax)

            plt.tight_layout()
            plt.savefig("{}{:05}.png".format(path,t),pad_inches=0)
            plt.clf() 

    logger.info("Movie...")
    if not os.path.exists("{0}video.mp4".format(path)):
        subprocess.call("avconv -y -r 3 -i {0}%05d.png {0}video.mp4".format(path), shell=True)
    subprocess.call("mplayer {}video.mp4  -idle -fixed-vo".format(path), shell=True)
