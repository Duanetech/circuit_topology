import numpy as np 

def export_psc(psclist):
    arr = np.array(psclist,dtype='object')
    np.savetxt('results/statistics/pscresults.csv',arr,fmt='%s,%s,%s,%s')