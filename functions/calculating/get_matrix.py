import numpy as np

def get_matrix(cmap3,protid):
    
    #Identify non zero values in the upper triangulated matrix of cmap3 and retrieve the respective indices
    rows,cols = np.nonzero(np.triu(cmap3))
    
    #create a numerical and character matrix based on the amount of nonzero values found in the previous function
    mat = np.identity(len(rows),dtype='int')*8
    c = np.empty([len(rows),len(rows)],dtype=object)

    #Change the values based on the type of connection
    for z in range(0,len(rows)):
        c[z][z] = '-'

    for x in range(0,len(rows)):
        for y in range(x+1,len(rows)):
            
            k = rows[x]
            m = cols[x]
            i = rows[y]
            j = cols[y]
            if k < i and m > j:
                mat[x][y] = 2
                mat[y][x] = 1
                c[x][y] = 'P-1'
                c[y][x] = 'P'
            elif m < i:
                mat[x][y] = 3
                mat[y][x] = 3
                c[x][y] = 'S'
                c[y][x] = 'S'
            elif k == i and m < j:
                mat[x][y] = 5
                mat[y][x] = 6
                c[x][y] = 'CP'
                c[y][x] = 'CP-1'
            elif k == i and m > j or m == j:
                mat[x][y] = 6
                mat[y][x] = 5
                c[x][y] = 'CP-1'
                c[y][x] = 'CP'
            elif m == i:
                mat[x][y] = 7
                mat[y][x] = 7
                c[x][y] = 'CS'
                c[y][x] = 'CS'
            elif m > i and j > m:
                mat[x][y] = 4
                mat[y][x] = 4
                c[x][y] = 'X'
                c[y][x] = 'X'
            else:
                print('Something went wrong')
            
    return mat, c