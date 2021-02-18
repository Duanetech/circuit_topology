import numpy as np


def secondary_struc_filter(cmap3,struc,filtered_structures = ['H','G'],ss_elements = ['H','E','B','G']):
    
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
            if i == struc_length:
                nstruc = nstruc + 1
            elif struc[i+1] != struc[i]:
                nstruc = nstruc + 1
    nstruc = nstruc - 1
    
    cmap4 = np.array(cmap3,copy=True)
    for i in range(len(cmap3)):
        for j in range(len(cmap3)):
            if cmap3[i][j] == 1:
                if struc_id[i] == struc_id[j] and struc[i] in filtered_structures:
                    cmap4[i][j] = 0

    return cmap4