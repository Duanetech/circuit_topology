import numpy as np

def energy_cmap(cmap3,numbering,res_names,potential_sign = 1):
    
    n_amino = 20
    aamap = np.array(['CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR','ALA','GLY','THR','SER','GLN','ASN','GLU','ASP','HIS','ARG','LYS','PRO'])
    y,x = np.nonzero(cmap3)
    
    numbering = np.array([numbering]).T

    res_names = np.array(res_names)
    resnames1 = res_names[y]
    resnames2 = res_names[x]

    potential_matrix = np.loadtxt('input_files/matrix_potential.txt','object')
    nan_matrix = potential_matrix == 'n'
    nan_matrix = nan_matrix * 1
    potential_matrix[potential_matrix == 'n'] = 100
    potential_matrix = potential_matrix.astype('float')

    for i in range(0,n_amino):
        for j in range(0,n_amino):
            if potential_matrix[i][j] == 0:
                potential_matrix[i][j] = 10
                potential_matrix[j][i] = 10
            if nan_matrix[i][j] == 1:
                potential_matrix[i][j]=potential_matrix[j][i]

    energy_cmap = np.zeros([len(cmap3),len(cmap3)],'float')

    for i in range(0,len(x)):

            energy_cmap[x[i]][y[i]] = potential_matrix[np.where(aamap==resnames1[i])[0][0]][np.where(aamap==resnames2[i])[0][0]]
            energy_cmap[y[i]][x[i]] = potential_matrix[np.where(aamap==resnames1[i])[0][0]][np.where(aamap==resnames2[i])[0][0]]
  
    if potential_sign == 1:
        energy_cmap = (energy_cmap > 0) * 1
    else:
        energy_cmap = (energy_cmap < 0) * 1
 
    return energy_cmap

