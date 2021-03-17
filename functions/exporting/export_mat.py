import numpy as np
import pandas as pd

def export_mat(cmap3,numbering,c,protid):
    numbering = np.array(numbering)
    rows,cols = np.nonzero(np.triu(cmap3))
    rowsnew = numbering[rows]
    colsnew = numbering[cols]
    names = []
    for i in range(len(rowsnew)):
        names.append(f'{rowsnew[i]} - {colsnew[i]}')
    df = pd.DataFrame(c,columns=names,index=names)
    df.to_csv('results/matrix/'+protid+'_mat.csv',sep=';')
