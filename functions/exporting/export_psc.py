"""
Created on Mon May 24 17:00:09 2021

@author: DuaneM

For exporting amount of PSC contacts to a csv file. 
"""
import numpy as np 

def export_psc(psclist):
    if len(psclist[0]) == 4:
        arr = np.array(psclist,dtype='object')
        np.savetxt('results/statistics/pscresults.csv',arr,fmt='%s,%s,%s,%s')
        print(f'Succesfully saved pscresults.csv')
        
    if len(psclist[0]) > 4:
        arr = np.array(psclist,dtype='object')
        np.savetxt('results/statistics/pscresults.csv',arr,fmt='%s,%s,%s,%s,%s,%s,%s')
        print(f'Succesfully saved pscresults.csv')