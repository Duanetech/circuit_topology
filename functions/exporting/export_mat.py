import numpy as np
import pandas as pd

def export_mat(index,mat,protid):
    
    d = dict({0:'-',1:'S',2:'P',3:'P-1',4:'C',5:'CP',6:'CP-1',7:'CS'})
    c = np.vectorize(d.get)(mat.astype(int))

    names = []
    for i in range(len(index)):
        names.append(f'[{index[i,0]} - {index[i,1]}]')
        
    df = pd.DataFrame(c,columns=names,index=names)
    df.to_csv('results/matrix/'+protid+'_mat.csv',sep=';')
