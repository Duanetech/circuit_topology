from  Bio.PDB.DSSP import DSSP
from Bio.PDB import Selection

def retrieve_secondary_struc(chain,input_path):
    
    # Uses biopython's built-in DSSP Function to retrieve the secondary structure
    # note! STRIDE and DSSP agree in 95,4% of the cases. DSSP tends to assign shorter secondary structures 
    # Must have DSSP installed --> see tutorial for instructions
    # conda install -c salilab dssp
    # https://en.wikipedia.org/wiki/STRIDE

    # H - Alpha-Helix
    # B - Isolated Beta-Bridge
    # E - Strand
    # G - 3-10 Helix
    # I - Pi helix
    # T - Turn
    # S - Bend

    model = chain.get_parent()
    res_list = Selection.unfold_entities(chain,"R")
    chain_len = len(res_list)
    
    dssp = DSSP(model,input_path)
    dssplist = list(dssp)[:chain_len]
    
    seq = [row[1] for row in dssplist]
    seq = ''.join(seq)
    struc = [row[2] for row in dssplist]
    struc = ''.join(struc)
    struc = struc.replace('-',' ')

    if len(struc) > len(seq):
        struc = struc[0:len(seq)]
    
    return seq,struc

