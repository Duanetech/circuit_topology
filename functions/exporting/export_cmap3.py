import numpy as np 

def export_cmap3(cmap3,protid):
    fname = 'results/circuit_diagram/'+protid+'-cmap3.csv'
    np.savetxt(fname,cmap3,fmt='%i',delimiter=',')