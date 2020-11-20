import numpy as np  

def get_stats(mat,protid):

    numparralel = (np.count_nonzero(mat == 1) + np.count_nonzero(mat == 2) + np.count_nonzero(mat == 5) + np.count_nonzero(mat == 6))/2
    numcross = np.count_nonzero(mat == 4)/2
    numseries = (np.count_nonzero(mat == 3) + np.count_nonzero(mat == 7))/2

    psc = np.array([numparralel,numseries,numcross],dtype='int')

    entangled = np.zeros([len(mat),1])

    for i in range(0,len(mat)-1):
        diag = np.diag(mat,k=i)
        entangled[i] = 1 - (sum(diag == 3)+ sum(diag == 7))/len(diag)
    
    return psc, entangled