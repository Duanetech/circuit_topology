from Bio.PDB import MMCIFParser, PDBParser
from Bio import BiopythonWarning
import warnings

def retrieve_chain(input_file):
    if input_file.endswith('cif'):
        
        input_file= 'input_files/mmcif/' + input_file
        #Supress harmless warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            #Import the protein data
            structure = MMCIFParser().get_structure(input_file.replace('.cif','').replace('input_files/mmcif/',''),input_file)
   
    else:
        input_file= 'input_files/pdb/' + input_file
        #Supress harmless warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            #Import the protein data
            structure = PDBParser(PERMISSIVE=1).get_structure(input_file.replace('.pdb','').replace('input_files/pdb/',''),input_file)
     
    chain = structure[0]['A']

    #remove heteroatoms (water molecules)
    removeres = []
    for res in chain:
        if res.id[0] != ' ':
            removeres.append(res.id)
            
    for res in removeres:
        chain.detach_child(res)
    return chain


