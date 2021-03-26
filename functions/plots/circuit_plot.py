import matplotlib.pyplot as plt
import numpy as np

def circuit_plot(index,protid,numbering):
    plt.ion()
    nseg = len(numbering)
    segment = list(range(0,nseg))  
    ax = plt.subplots()[1]
    plt.plot()
    plt.xlim(0,len(numbering)+1)
    plt.hlines(0,0,len(numbering)+1,colors = 'k',linestyles='dashed')
    ax.tick_params(labelleft = False)
    ax.set_xlabel('Residue')

    for row in index:
        X = np.array([row[0],row[1]])
        r = (X[1]-X[0])/2
        center = np.mean(X)
        xx = np.arange(X[0],X[1]+0.001,0.01)
        np.seterr(invalid = 'ignore')
        yy = np.sqrt(r**2 - (xx - center)**2)
        where_nans = np.isnan(yy)
        yy[where_nans] = 0
        plt.plot(xx,yy,
                color ='k')