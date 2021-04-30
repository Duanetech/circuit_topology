from Bio.PDB import MMCIFParser, PDBParser
from Bio import BiopythonWarning
import warnings

def retrieve_chain(input_file,chainid = 0):
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

    model = structure[0]
    chainlist = model.get_list()

    residue_to_remove = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != ' ':
                residue_to_remove.append((chain.id, residue.id))


    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    if type(chainid) == int:
        chain = model.get_list()[chainid]
    elif type(chainid) == str:
        chain = model[chainid]
    else:
        raise TypeError

    protid = structure.id+ '_' + chain.id

    return chain,protid


