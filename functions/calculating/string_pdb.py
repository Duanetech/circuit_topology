import numpy as np
import copy

def string_pdb(cmap,sitelist,numbering,threshold):
    #% this gives a reduced form of S and d, where nodes containing no contacts are excluded
    cmap = (cmap>0) * 1
    
    for i in range(0,len(cmap)):
        cmap[i][i] = 0

    y,x = np.nonzero(cmap)

    S = np.zeros([len(y),len(y)],dtype='int')
    cnt = np.ones([1,len(cmap)],dtype='int')


    for i in range(0,len(y)):
        a1 = np.array(np.where(y==x[i]))
        n = cnt[0][x[i]]
        cnt[0][x[i]] = cnt[0][x[i]] + 1
        S[i][a1[0][n-1]] = 1
    d = y
    # simplify matrix (remove 0-0 lines and rows)
    y,x = np.nonzero(S)
    
    #transform contact indices in one minima vector
    a1 = np.minimum(x,y)
    s1 = copy.deepcopy(a1)
    s2 = copy.deepcopy(a1) 
    
    #produce strings from the minima vector
    i = 0

    while s2.size != 0:
        s1[np.where(s1==np.amin(s2))] = i
        s2 = np.setdiff1d(s2,np.amin(s2))
        i = i + 1

    b1 = (a1 == x)
    parentheses = np.zeros(len(d),dtype='int')
    j = 1

    for i in range(0,max(d)+1):

        c1 = d[d==i]
        c2 = b1[d==i]
        c3 = s1[d==i]

        if len(c1) > 1:
            s1[d==i] = np.concatenate([np.flipud(c3[c2]),np.flipud(c3[~c2])])
            parentheses[d==i] = j
            j = j + 1

    #delete long-range contacts

    s2 = s1
    dell = np.zeros(len(s1),dtype='int')

    for i in range(0,len(s1)):
        if (sitelist[s1[i]][1]-sitelist[s1[i]][0]) > threshold:
            dell[i] = 1
    s2 = np.delete(s2,np.where(dell != 0))

    #look for circuits
    n = np.zeros([len(s2)],dtype='int')
    cnt = np.zeros([1,max(s2)+1],dtype='int')
    j = 1


    for i in range(0, len(s2)):
        n[i] = j
        cnt[0][s2[i]] = cnt[0][s2[i]] + 1
        if np.mean(cnt[cnt>0]) == 2:
            j = j + 1
            cnt = np.zeros([1,max(s2)+1],dtype='int')
        i = i + 1

    #combine consecutive contacts (e.g., 223344)
    num = []
    for i in range(1,np.amax(n)+1):
        num.append(np.count_nonzero(n==i))

    c = np.zeros(sum(num),dtype='int')

    for i in range(0,len(num)):
        for j in range(0,len(c)):
            if num[i] > 2:
                c[j] = 0
            else:
                c[j] = 1

    d = n * c

    for i in range(0,len(d)):
        if d[i] > 0:
            if i > 0 and d[i-1] > 0:
                d[i] = d[i-1]

    n[d>0] = d[d>0]      
    #n2 = np.unique(n)

    segends = []
    for i in range(1,max(n)+1):
        segends.append(sitelist[s2[max(np.where(n==i)[0])]][1])

    segnums = max(n)
    length = np.array(segends[0])
    length = np.append(length,np.diff(segends)[0])
    meanlength = np.mean(length)
    standard_dev = np.std(length,ddof=1)
    threshold_length = meanlength - standard_dev / 2

    #delete those circuits that are too small
    index_circuit = np.zeros(segnums)
    for i in range(0,len(length)):
        if segnums != 1 and length[i] <= threshold_length:
            index_circuit[i] = 1

    segnums = int(segnums - np.sum(index_circuit))

    return segnums,meanlength,segends