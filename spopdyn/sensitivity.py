import numpy as np

import matplotlib.pyplot as plt


def sensitivity_analysis(D, X, Y=None):
    """
    Uses Taylor's theorem to estiamte the sensitivity of the system to
    parameters. 

    Y = X + DS <=> S = D^-1 (Y-X) = S
    """
    if Y is None:
        dX = X
    else:
        dX = Y-X

    # make a column vector if it is not the case.
    if len(dX.shape) == 1:
        dX = dX.reshape(-1,1)

    if len(D.shape) == 1:
        D = np.vstack([D**n for n in  range(len(D))])


    print D
    Dinv = np.linalg.inv(D)
    print dX, Dinv
    S = np.dot(Dinv,dX)
    
    return S

def plot_poly(S,x=None):
    if x is None:
        x = np.linspace(-1,1,200)
    y = np.poly1d(S[::-1])(x)
    plt.plot(x,y)
