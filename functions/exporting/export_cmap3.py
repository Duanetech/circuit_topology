import numpy as np 

def export_cmap3(index,protid,numbering):
    cmap = np.zeros([len(numbering),len(numbering)],dtype='int')
    
    for row in index:
        x = row[0]
        y = row[1]
        cmap[x][y] = 1
        cmap[y][x] = 1

    fname = 'results/circuit_diagram/'+protid+'-cmap3.csv'
    np.savetxt(fname,cmap,fmt='%i',delimiter=',')