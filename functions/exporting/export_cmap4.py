"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

For transforming Segment contact map indices to a contact map and exporting that to a csv file.
"""
import numpy as np
import pandas as pd 

def export_cmap4(cmap4,segment,structure,protid):
    names = []

    for i in range(1,max(segment)+1):
       idx = structure[np.where(segment == i)[0][0]]
       names.append(str(i)+' - '+idx)

    df = pd.DataFrame(cmap4,columns=names,index=names)
    df.to_csv('results/circuit_diagram/'+protid+'_cmap4.csv',sep=';')
    print(f'Succesfully saved {protid}_cmap4.csv')