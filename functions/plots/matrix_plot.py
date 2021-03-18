import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import warnings

def matrix_plot(mat,protid):

    #create custom colormap
    newcolors = np.array([[218/255, 219/255, 228/255,1],
                      [131/255, 139/255, 197/255,1],
                      [172/255,200/255,247/255,1], 
                      [174/255,213/255,129/255,1], 
                      [186/255, 155/255, 201/255,1], 
                      [172/255, 200/255, 247/255,1],
                      [174/255, 213/255, 129/255,1],
                      [131/255,139/255, 197/255,1]])
    newcmp = ListedColormap(newcolors)

    fig, ax = plt.subplots()
    color = plt.get_cmap(newcmp, 8)

    #plot data
    pngmat = plt.imshow(mat,cmap=color,vmin = np.min(mat)-.5, vmax = np.max(mat)+.5)
    ax.set_title(protid)
    ax.tick_params(labelleft = False,labelbottom = False,bottom = False,left= False)
    cbar = fig.colorbar(pngmat)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cbar.ax.set_yticklabels(['-','S','P','P-1','X','CP','CP-1','CS'])
