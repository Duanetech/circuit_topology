from Bio.PDB import Selection,NeighborSearch
from collections import Counter
import numpy as np

def get_cmap(
            chain,
            cutoff_distance = 4.5,
            cutoff_numcontacts = 5,        
            exclude_neighbour = 3):


    atom_list = Selection.unfold_entities(chain,'A')
    res_list = Selection.unfold_entities(chain,'R')

    ns = NeighborSearch(atom_list)
    all_neighbours = ns.search_all(cutoff_distance,'A')

    res_names = []
    numbering = []
    
    for res in res_list:
        res_names.append(res.get_resname())
        numbering.append(res.get_id()[1])

    numbering = np.array(numbering)
    segment = np.array(range(len(numbering)))

    index_list = []
    for atompair in all_neighbours:
        res1 = segment[numbering == atompair[0].get_parent().id[1]][0]
        res2 = segment[numbering == atompair[1].get_parent().id[1]][0]

        if abs(res1-res2) > exclude_neighbour:
            index_list.append((res1,res2))

    index_list.sort()
    count = Counter(index_list)

    new_index = [values for values in count if count[values] >= cutoff_numcontacts]
    
    return np.array(new_index),numbering,res_names