import numpy as np  

def get_stats(mat):
    entangled = np.zeros([len(mat),1])

    for i in range(0,len(mat)-1):
        diag = np.diag(mat,k=i)
        entangled[i] = 1 - (sum(diag == 1)+ sum(diag == 7))/len(diag)
    
    return entangled