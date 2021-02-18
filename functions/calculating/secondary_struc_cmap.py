import numpy as np
from Bio.PDB import PDBParser,Selection
from scipy.spatial.distance import pdist, squareform

def secondary_struc_cmap(chain,seq,struc,cutoff_distance = 6,cutoff_numcontacts = 10,
            exclude_neighbour=0,ss_elements = ['H','E','B','G']):

    struc_length = len(struc)       
    
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

    for res in res_list:
        res_names.append(res.get_resname())

    for num, atom in enumerate(atom_list):
        residue_number[num] = atom.get_parent().get_id()[1]
        coords[num] = atom.get_coord()
        atom_names.append(atom.get_name())
        y[num] = num

    numbering = list(range(residue_number[0],residue_number[-1]+1))
    residue_number = residue_number - numbering[0] 
    natoms = len(residue_number)       

    nseg = 1
    segment = np.zeros([struc_length],dtype='int')

    for i in range(0,struc_length):
        if struc[i] in ss_elements:
            segment[i] = nseg
            if i == struc_length:
                nseg = nseg + 1
            elif struc[i+1] != struc[i]:
                nseg = nseg + 1
    nseg = nseg - 1

    #atom-atom based contact map
    cmap = squareform(pdist(coords))
    cmap = (cmap < cutoff_distance) * 1

    #Secondary structure based contact map
    cmap2 = np.zeros([nseg,nseg],dtype='int')

    for i in range(0,natoms):
        for j in range(i+1,natoms):
            diff_res = abs(residue_number[i]-residue_number[j])
            
            if cmap[i][j] == 1  and diff_res != 1 and diff_res != 2 and diff_res != 3:
                seg_i = segment[residue_number[i]]
                seg_j = segment[residue_number[j]]
            
                if seg_i != 0 and seg_j != 0 :
                    cmap2[seg_i-1][seg_j-1] = cmap2[seg_i-1][seg_j-1]+1
                        
    cmap2 = cmap2 + cmap2.T

    #set values close to diagonal to zero
    for i in range(0,exclude_neighbour+1):
        for j in range(0,len(cmap2)-i):
            cmap2[j][j+i] = 0
            cmap2[j+i][j] = 0

    cmap5 = (cmap2 >= cutoff_numcontacts) * 1
    
    return cmap5,cmap,cmap2,struc,segment,numbering