import numpy as np 
import pandas as pd 

def export_cmap3(index,protid,numbering):
    cmap = np.zeros([len(numbering),len(numbering)],dtype='int')
    
    for row in index:
        x = row[0]
        y = row[1]
        cmap[x][y] = 1
        cmap[y][x] = 1

    df = pd.DataFrame(cmap,index=numbering,columns=numbering)
    df.to_csv('results/circuit_diagram/'+protid+'-cmap3.csv') 