import numpy as np
from Bio.PDB import PDBParser,Selection
from scipy.spatial.distance import pdist, squareform

def secondary_struc_cmap(chain,seq,struc,cutoff_distance = 4.5,cutoff_numcontacts = 10,
            exclude_neighbour=0,ss_elements = ['H','E','B','G']):

    struc_length = len(struc)       
    
    protid = chain.get_parent().get_parent().id

    #Get a list of the residues and atoms
    res_list = Selection.unfold_entities(chain,"R")
    atom_list = Selection.unfold_entities(chain,"A")

    #Make list of the atom information
    
    #Residue information
    res_names = []
    numbering = []

    #Atom information
    residue_number = np.zeros(len(atom_list),dtype='int')
    atom_names = []
    coords = np.zeros([len(atom_list),3])

    for num,res in enumerate(res_list):
        res_names.append(res.get_resname())
        numbering.append(res.get_id()[1])

        
    for num, atom in enumerate(atom_list):
        residue_number[num] = atom.get_parent().get_id()[1]
        coords[num] = atom.get_coord()
        atom_names.append(atom.get_name()) 

    residue_number_new = residue_number - numbering[0]
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

    #Delete duplicates
    duplicate = np.zeros(len(atom_list),dtype='int')
    for i in range(1,len(duplicate)):
        if residue_number[i] == residue_number[i-1] and atom_names[i] == atom_names[i-1]:
            duplicate[i] = 1
    
    residue_number = residue_number[np.where(duplicate != 1)]
    coords = coords[np.where(duplicate != 1)]
    natoms = len(residue_number)
    
    #atom-atom based contact map
    cmap = (squareform(pdist(coords)) < cutoff_distance) * 1
    np.fill_diagonal(cmap,0)
    assert len(seq) == len(numbering),f'PDB file and Secondary structure map do not match!\n {protid} - PDB: {len(numbering)} Residues VS. SS: {len(seq)} Residues. '

    #Secondary structure based contact map
    cmap2 = np.zeros([nseg,nseg],dtype='int')

    for i in range(0,natoms):
        for j in range(i+1,natoms):
            diff_res = abs(residue_number[i]-residue_number[j])
            
            if cmap[i][j] == 1  and diff_res != 1 and diff_res != 2 and diff_res != 3:
                seg_i = segment[residue_number_new[i]]
                seg_j = segment[residue_number_new[j]]
            
                if seg_i != 0 and seg_j != 0 :
                    cmap2[seg_i-1][seg_j-1] = cmap2[seg_i-1][seg_j-1]+1
                        
    cmap2 = cmap2 + cmap2.T

    #set values close to diagonal to zero
    for i in range(0,exclude_neighbour+1):
        for j in range(0,len(cmap2)-i):
            cmap2[j][j+i] = 0
            cmap2[j+i][j] = 0

    cmap4 = (cmap2 >= cutoff_numcontacts) * 1
    
    return cmap4,segment,numbering,protid