import numpy as np
import pandas as pd 

def export_cmap4(cmap4,segment,structure,protid):
    names = []

    for i in range(1,max(segment)+1):
       idx = structure[np.where(segment == i)[0][0]]
       names.append(str(i)+' - '+idx)

    df = pd.DataFrame(cmap4,columns=names,index=names)
    df.to_csv('results/circuit_diagram/'+protid+'_cmap4.csv',sep=';')