from Bio.PDB import MMCIFParser, PDBParser
from Bio import BiopythonWarning
import warnings

def retrieve_chain(input_file):
    if input_file.endswith('cif'):
        
        input_filepath= 'input_files/cif/' + input_file
        #Supress harmless warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            #Import the protein data
            structure = MMCIFParser().get_structure(input_file.replace('.cif',''),input_filepath)
   
    else:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
        
            input_filepath= 'input_files/pdb/' + input_file

            structure = PDBParser(PERMISSIVE=1).get_structure(input_file.replace('.pdb',''),input_filepath)

    chainlist = structure[0].get_list()
    chain = chainlist[0]
    
    #remove heteroatoms (water molecules)
    removeres = []
    for res in chain:
        if res.id[0] != ' ':
            removeres.append(res.id)
            
    for res in removeres:
        chain.detach_child(res)
        
    return chain


