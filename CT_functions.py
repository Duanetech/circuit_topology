from Bio import BiopythonWarning
from Bio.PDB import *
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.colors import ListedColormap
from scipy.spatial.distance import pdist, squareform
import sys
import string
import warnings
import os

def retrieve_pdb_files():
    server = PDBList(server='ftp://ftp.wwpdb.org', pdb=None, obsolete_pdb=None ,verbose=True)
    pdb_list = open('PDB/pdblist.txt','r')
    content = pdb_list.read().split()
    pdb_list.close()
    server.download_pdb_files(content,pdir="PDB/mmCif",file_format='mmCif', overwrite=True,obsolete= False)

def circuit_diagram_residue(cif_file,cutoff_distance = 3.6,cutoff_numcontacts = 3,plot_figs = 1,fileformat='jpeg'):
    
    #Supress harmless warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        #Import the protein data
        
        structure = MMCIFParser().get_structure(cif_file,cif_file)
        chain = structure[0]['A']
        protid1 = cif_file.replace('.cif','')
        protid = protid1.replace('PDB/mmCif/','')
    
    #remove heteroatoms (water molecules)
    removeres = []
    for res in chain:
        if res.id[0] != ' ':
            removeres.append(res.id)
    for res in removeres:
        chain.detach_child(res) 
    
    #Get a list of the residues and atoms
    res_list = Selection.unfold_entities(chain,"R")
    atom_list = Selection.unfold_entities(chain,"A")

    #Make list of the atom information
    residue_number = np.zeros(len(atom_list),dtype='int')
    name = []
    coords = np.zeros([len(atom_list),3])
    y = np.zeros(len(atom_list),dtype='int')

    for num, atom in enumerate(atom_list):
        residue_number[num] = atom.get_parent().get_id()[1]
        coords[num] = atom.get_coord()
        name.append(atom.get_name())
        y[num] = num
    numbering = list(range(residue_number[0],residue_number[-1]+1))
    residue_number = residue_number - numbering[0] 

    #divide into segments based on residue
    nseg = len(res_list)
    segment = list(range(0,nseg))

    #delete duplicate atoms due to multiple occupancy
    duplicate = np.zeros(len(atom_list),dtype='int')
    for i in range(1,len(duplicate)):
        if residue_number[i] == residue_number[i-1] and name[i] == name[i-1]:
            duplicate[i] = 1


    residue_number = residue_number[np.where(duplicate != 1)]
    coords = coords[np.where(duplicate != 1)]

    #atom-atom based contact map
    cmap = squareform(pdist(coords))
    natoms = len(atom_list)

    cmap[cmap < cutoff_distance]= 1
    cmap[cmap > cutoff_distance]= 0

    #segment-segment based contact map, based on atom-atom based contact map
    cmap2 = np.zeros([nseg,nseg],dtype='int')

    for i in range(0,natoms):
        for j in range(i+1,natoms):
            
            seg_i = segment[residue_number[i]]
            seg_j = segment[residue_number[j]]
            if cmap[i][j] == 1:
                if seg_i != 0 and seg_j != 0:
                    cmap2[seg_i][seg_j] = cmap2[seg_i][seg_j]+1

    cmap2 = cmap2 + cmap2.T

    #set values close to diagonal to zero
    for i in range(0,len(cmap2)):
        cmap2[i][i] = 0
    for i in range(0,len(cmap2)-1):
        cmap2[i][i+1] = 0
        cmap2[i+1][i] = 0
    for i in range(0,len(cmap2)-2):
        cmap2[i][i+2] = 0
        cmap2[i+2][i] = 0
    for i in range(0,len(cmap2)-3):
        cmap2[i][i+3] = 0
        cmap2[i+3][i] = 0
    #plot figures
    if plot_figs == 1:
        fig, ax = plt.subplots()
        plt.plot()
        plt.hlines(0,0,len(numbering)+1,colors = 'k')
        plt.xlim(1,len(numbering)+1)
        ax.tick_params(labelleft = False)
        ax.set_xlabel('Residue')
        ax.set_title(protid)
        for i in range(0,nseg):
            for j in range(i + 1,nseg):
                if cmap2[i][j] >= 3:
                    X = np.array([segment[i]+1,segment[j]+1])
                    r = (X[1]-X[0])/2
                    center = np.mean(X)
                    xx = np.arange(X[0],X[1]+0.001,0.01)
                    np.seterr(invalid = 'ignore')
                    yy = np.sqrt(r**2 - (xx - center)**2)
                    where_nans = np.isnan(yy)
                    yy[where_nans] = 0
                    plt.plot(xx,yy,
                    color ='k',
                    linewidth = np.log(cmap2[i][j]+1))
        plt.savefig('PDB/Contact map/' + protid + '_cmap.'+ fileformat)

    #transform the residue based contact map based on the minimal number of contacts
    cmap3 = (cmap2 >= cutoff_numcontacts) * 1
    
    return cmap3, protid , cmap, numbering 

def get_matrix(cmap3,protid,plot_figs=1,fileformat='jpeg'):
    
    #Identify non zero values in the upper triangulated matrix of cmap3 and retrieve the respective indices
    rows,cols = np.nonzero(np.triu(cmap3))
    
    #create a numerical and character matrix based on the amount of nonzero values found in the previous function
    mat = np.identity(len(rows))*8
    c = np.empty([len(rows),len(rows)],dtype=object)

    #Change the values based on the type of connection
    for z in range(0,len(rows)):
        c[z][z] = 'P'

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
    if plot_figs:
        #create custom colormap
        newcolors = np.array([[172/255,200/255,247/255,1], [174/255,213/255,129/255,1], [131/255,139/255, 197/255,1], [186/255, 155/255, 201/255,1], [172/255, 200/255, 247/255,1], [174/255, 213/255, 129/255,1], [131/255, 139/255, 197/255,1], [218/255, 219/255, 228/255,1]])
        newcmp = ListedColormap(newcolors)

        fig, ax = plt.subplots()
        color = plt.get_cmap(newcmp, 8)

        #plot data
        pngmat = plt.imshow(mat,cmap=color,vmin = np.min(mat)-.5, vmax = np.max(mat)+.5)
        ax.set_title(protid)
        plt.clim(0.5,8.5)
        ax.tick_params(labelleft = False,labelbottom = False,bottom = False,left= False)
        cbar = fig.colorbar(pngmat)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cbar.ax.set_yticklabels(['P','P-1','S','X','CP','CP-1','CS','-'])
        plt.savefig('PDB/Matrix/' +protid + '_matrix.'+ fileformat)
    return mat, c

def get_ct_stats(mat,protid,plot_figs=1,fileformat='jpeg'):

    numparralel = (np.count_nonzero(mat == 1) + np.count_nonzero(mat == 2) + np.count_nonzero(mat == 5) + np.count_nonzero(mat == 6))/2
    numcross = np.count_nonzero(mat == 4)/2
    numseries = (np.count_nonzero(mat == 3) + np.count_nonzero(mat == 7))/2

    psc = np.array([numparralel,numseries,numcross],dtype='int')

    entangled = np.zeros([len(mat),1])

    for i in range(0,len(mat)-1):
        diag = np.diag(mat,k=i)
        entangled[i] = 1 - (sum(diag == 3)+ sum(diag == 7))/len(diag)
    
    if plot_figs:
        fig,(ax1,ax2)= plt.subplots(1,2)
        ax1.plot(entangled)
        ax2.pie(psc,autopct = '%1.1f%%')
        ax1.set_xlabel('Distance from diagonal')
        ax1.set_ylabel('Fraction entangled')
        fig.suptitle(protid)
        ax2.legend(['Paralell','Series','Cross'],bbox_to_anchor=(.5, -0.5, 0.5, 0.5) )
        plt.savefig('PDB/Statistics/' +protid + '_ctstats.'+fileformat)

    return psc, entangled




