'''
from CT_functions import retrieve_pdb_files, circuit_diagram_residue, get_matrix, get_ct_stats
import os

#Enter a list with protein names that need to be downloaded (optional)
dlinput = input('Do you want to fetch PDB\'s from Database? (y/n) ')
if dlinput == 'y':  
    retrieve_pdb_files()



for files in os.listdir('input_files/mmcif'):
    
    cif_file = 'input_files/mmcif/'+ files
    
    #Step 1 - Draw a segment-segment based contact map
    cmap3, protid, cmap, numbering = circuit_diagram_residue(cif_file)
 
    #Step 2 - Draw a circuit topology relations matrix
    mat, c = get_matrix(cmap3,protid)

    #Step 3 - Circuit topology statistics
    psc, entangled = get_ct_stats(mat,protid)





'''