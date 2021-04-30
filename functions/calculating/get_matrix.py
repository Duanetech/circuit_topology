import numpy as np

def get_matrix(index):
    
    #create a numerical and character matrix based on the amount of nonzero values found in the previous function
    mat = np.zeros((len(index), len(index)),dtype = 'int')

    #Change the values based on the type of connection
    

    P = 0
    S = 0
    X = 0

    for x in range(0,len(index)):
        i = index[x,0]
        j = index[x,1]
        for y in range(x+1,len(index)):
            k = index[y,0]
            l = index[y,1]
            #series
            if (j < k):
                S=S+1
                mat[x, y]=1
                mat[y, x]=1

            #parallel    
            elif (i>k and j<l):
                P=P+1
                mat[x, y]=2
                mat[y, x]=3
            
            #5: CP
            #6: CP-1    
            elif (i==k and j<l):
                mat[x, y]=5
                mat[y, x]=6
                P += 1
           
            elif (i==k and l<j):
                mat[x, y]=6
                mat[y, x]=5
                P += 1
               
            elif (k>i and j==l):
                mat[x,y]=6
                mat[y,x]=5
                P += 1
                
            elif(i>k and l==j):
                mat[x,y]=5
                mat[y,x]=6
                P += 1
            #inverse parallel
            elif (k>i and l<j):
                P += 1
                mat[x, y]=3
                mat[y, x]=2
            #CS
            elif (j ==k):
                mat[x,y]=7
                mat[y,x]=7
                S += 1
            
            if (k>i and k<j and j<l):
                X += 1
                mat[x, y]=4
                mat[y, x]=4
            elif (i>k and i< l and j> l):
                X += 1
                mat[x, y]=4
                mat[y, x]=4

    psc = np.array([P,S,X],dtype='int')

    return mat,psc