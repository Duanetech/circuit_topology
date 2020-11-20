import matplotlib.pyplot as plt
import numpy as np

def circuit_plot(cmap,protid,numbering,fileformat='jpeg'):
    nseg = len(cmap)
    segment = list(range(0,nseg))  
    ax = plt.subplots()[1]
    plt.plot()
    plt.hlines(0,0,len(numbering)+1,colors = 'k')
    plt.xlim(1,len(numbering)+1)
    ax.tick_params(labelleft = False)
    ax.set_xlabel('Residue')
    ax.set_title(protid)
    for i in range(0,nseg):
        for j in range(i + 1,nseg):
            if cmap[i][j] >= 3:
                X = np.array([segment[i]+1,segment[j]+1])
                r = (X[1]-X[0])/2
                center = np.mean(X)
                xx = np.arange(X[0],X[1]+0.001,0.01)
                np.seterr(invalid = 'ignore')
                yy = np.sqrt(r**2 - (xx - center)**2)
                where_nans = np.isnan(yy)
                yy[where_nans] = 0
                plt.plot(xx,yy,
                color ='k',
                linewidth = np.log(cmap[i][j]+1))
    plt.savefig('results/circuit_diagram/' + protid + '_cmap.'+ fileformat)
