import numpy as np


def secondary_struc_filter(index,struc,filtered_structures = ['H','G'],ss_elements = ['H','E','B','G','T','S','I']):
    
    #STRIDE
    # H - Alpha-Helix
    # E - Extended Configuration (Beta-sheet)
    # B - Isolated Beta Bridge
    # b - Isolated Beta Bridge
    # G - 3-10 Helix

    #DSSP
    # H - Alpha-Helix
    # B - Isolated Beta-Bridge
    # E - Strand
    # G - 3-10 Helix
    # I - Pi helix
    # T - Turn
    # S - Bend

    struc_length = len(struc)

    nstruc = 1
    struc_id = np.zeros([struc_length],dtype='int')

    for i in range(0,struc_length):
        if struc[i] in ss_elements:
            struc_id[i] = nstruc
            if i == struc_length-1:
                nstruc = nstruc + 1
            elif struc[i+1] != struc[i]:
                nstruc = nstruc + 1
    nstruc = nstruc - 1
    
    remove = []

    for num,i in enumerate(index):
        res1 = i[0]
        res2 = i[1]
        if struc_id[res1] == struc_id[res2] and struc[res1] in filtered_structures:
            remove.append(num)

    index_filtered = np.delete(index,remove,0)

    return index_filtered,nstruc