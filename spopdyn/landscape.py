from __future__ import division
import numpy as np


def seq_gaussian(size, alpha,rmax=1,compute_semivariogram=False):
    """
    Create a random landscape using the sequential gaussian algorithm.
    
    Args:
        size (int): number of gridpoints on one axis.
        alpha (float): spatial correlation
    Returns:
        (np.array) Landscape. 
    """
    N = size**2
    landscape = np.zeros((size,size))

    
    path = np.arange(N)
    np.random.shuffle(path)

    landscape.flat[path[0]] = np.random.uniform() 

    #print "Compute Torus distance..."
    D = np.zeros((N,N))
    for i,x in enumerate(path):
        for j,y in enumerate(path):
            D[i,j] = torus_distance(x,y,size)

    #print "Sequential gaussian"        
    #print "D: {}".format(D.shape)
    for i,p in enumerate(path[1:]):

        #print "{}th visited point: {}".format(i,p)
        K = landscape.flat[path[:i+1]]
        d = D[:i+2,:i+2]
        mu,sigma = gaussian_kriging_step(K,d,alpha,rmax)
        #print("{:0.2%} p:{}, mu: {}, sigma: {}".format(i/float(N),p,mu,sigma))
        landscape.flat[p] = np.random.normal(mu,sigma)

    
    
    if compute_semivariogram:
        semivariogram = np.zeros(D.max()+1)
        semi_n= np.zeros(D.max()+1)
        for i,x in enumerate(path):
            for j,y in enumerate(path):
                semivariogram[D[i,j]] += (landscape.flat[x]-landscape.flat[y])**2
                semi_n[D[i,j]] += 1
        semivariogram /= 2*semi_n
        return landscape,semivariogram
    else:
        return landscape

def torus_distance(a,b,size):
    """
    return the torus distance given the position in the flat_array
    """
    dx = abs(int(a/size) - int(b/size))
    dy = abs(a%size - b%size)

    # Torus !
    dx = min(dx,size-dx)
    dy = min(dy,size-dy)

    return np.sqrt(dx**2+dy**2)
    
def gaussian_correlogram(d,alpha,rmax=1):
    return rmax*np.exp(-3*d**2*alpha**-2)
def gaussian_semivariogram(d,alpha,rmax=1):
    return rmax*(1-np.exp((-3*d**2)/alpha**2))


def gaussian_kriging_step(Z,d,alpha,rmax):
    """
    Given the k first points, returns the distribution of the k+1 point.
    Args:
        Z (array): The value in the k first points
        d (matrix): k+1*k+1 distance matrix
        alpha (float): spatial autocorellation (for the semi-variogram)
        rmax (float): maximum correlation (for the semi-variogram)
    
    Returns:
        (tuple of float): mean and variance of the estimate of Z_k+1. 
    """
    k = len(Z)
    
    R = gaussian_correlogram(d[:-1,:-1],alpha,rmax)
    Rkp1 = gaussian_correlogram(d[-1,:-1],alpha,rmax)
    #print "d: {}, R: {}, Rkp1: {}".format(d.shape,R.shape,Rkp1.shape)
    
    W = np.dot(np.linalg.inv(R),np.transpose(Rkp1))
    
    mu = np.dot(Z,W)
    sigma = 1 - np.dot(W,Rkp1)

    return mu,sigma


def gradient(size,sigma=0.1):
    lines = []
    for i in range(size):
        lines.append(np.random.normal(i/size,sigma,size))
    return np.array(lines)

def main():
    import matplotlib.pyplot as plt
    size = 25
    plt.imshow(gradient(size),aspect="equal",interpolation="none")
    plt.show()
    
def test_gaussian():
    import matplotlib.pyplot as plt
    alpha = 3
    rmax = 1
    size = 10
    ld,sm = seq_gaussian(size=size,alpha=alpha,rmax=rmax)
    
    x = np.linspace(0,size)
    plt.subplot(2,1,1)
    plt.plot(x,gaussian_semivariogram(x,alpha,rmax=rmax),label="Semi variogramme Gaussien")
    plt.plot(range(len(sm)),sm,".-",label="Semi variogramme experimental")
    plt.legend(loc="lower right")
    plt.subplot(2,1,2)
    plt.imshow(ld,aspect="equal",interpolation="none")
    plt.show()

#main()










