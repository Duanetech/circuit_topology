from Bio.PDB import MMCIFParser,Selection
from Bio import BiopythonWarning
from scipy.spatial.distance import pdist, squareform
import numpy as np
import warnings

def get_cmap(chain,cutoff_distance = 3.6,cutoff_numcontacts = 3,length_filtering = 0):
    
    protid = chain.get_parent().get_parent().id
    
    #Get a list of the residues and atoms
    res_list = Selection.unfold_entities(chain,"R")
    atom_list = Selection.unfold_entities(chain,"A")
    
    #Make list of the atom information
    residue_number = np.zeros(len(atom_list),dtype='int')
    res_names = []
    atom_names = []
    coords = np.zeros([len(atom_list),3])
    y = np.zeros(len(atom_list),dtype='int')
    res_names = []

    for res in res_list:
        res_names.append(res.get_resname())

    for num, atom in enumerate(atom_list):
        residue_number[num] = atom.get_parent().get_id()[1]
        coords[num] = atom.get_coord()
        atom_names.append(atom.get_name())
        y[num] = num
        
    numbering = list(range(residue_number[0],residue_number[-1]+1))
    residue_number = residue_number - numbering[0] 
    
    #divide into segments based on residue
    nseg = len(numbering)
    segment = list(range(0,nseg+1))

    #delete duplicate atoms due to multiple occupancy
    duplicate = np.zeros(len(atom_list),dtype='int')
    for i in range(1,len(duplicate)):
        if residue_number[i] == residue_number[i-1] and atom_names[i] == atom_names[i-1]:
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
            if length_filtering != 0: 
                dist_sequence = abs(seg_i-seg_j)
                if cmap[i][j] == 1 and dist_sequence < length_filtering:
                    if seg_i != 0 and seg_j != 0:
                        cmap2[seg_i][seg_j] = cmap2[seg_i][seg_j]+1
            else:
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

    cmap3 = (cmap2 >= cutoff_numcontacts) * 1
    return cmap3, cmap2, protid ,numbering, res_names



