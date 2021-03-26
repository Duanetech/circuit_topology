from Bio.PDB import MMCIFParser,Selection
from Bio import BiopythonWarning
from scipy.spatial.distance import pdist, squareform
from collections import Counter
import numpy as np
import warnings

def get_cmap(chain,cutoff_distance = 4.5,cutoff_numcontacts = 5,
            length_filtering = 0,exclude_neighbour=3):

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

    for res in res_list:
        res_names.append(res.get_resname())
        numbering.append(res.get_id()[1])

    for num, atom in enumerate(atom_list):
        residue_number[num] = atom.get_parent().get_id()[1]
        coords[num] = atom.get_coord()
        atom_names.append(atom.get_name()) 
    
    #divide into segments based on residue
    nseg = len(res_list)
    segment = list(range(0,nseg+1))

    #delete duplicate atoms due to multiple occupancy
    duplicate = np.zeros(len(atom_list),dtype='int')
    for i in range(1,len(duplicate)):
        if residue_number[i] == residue_number[i-1] and atom_names[i] == atom_names[i-1]:
            duplicate[i] = 1
    
    residue_number = residue_number[np.where(duplicate != 1)]
    coords = coords[np.where(duplicate != 1)]

    #get indices of atom-atom contacts
    cmap = (squareform(pdist(coords)) < cutoff_distance) * 1
    np.fill_diagonal(cmap,0)
    atom_index = np.transpose(np.nonzero(np.triu(cmap)))

    #segment-segment based index, based on atom-atom based index

    res_index = []
    for i in atom_index:
        res_1 = np.where(numbering == residue_number[i[0]])[0][0]
        res_2 = np.where(numbering == residue_number[i[1]])[0][0]
        if abs(res_1-res_2) > exclude_neighbour:
            res_index.append((res_1,res_2))

    count = Counter(res_index)
    index = []
    for values in count:
        if count[values] >= 5:
            index.append(values)
   

    return np.array(index), numbering, protid ,res_names


