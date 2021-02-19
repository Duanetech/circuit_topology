from ipywidgets import widgets, HBox

def pdbbutton():
    fileformat = 'pdb'

def cifbutton():
    fileformat = 'cif'
    print('Do you want to fetch CIF\'s from Database? pdblist.txt')
    yes = widgets.Button(description='Yes')
    no = widgets.Button(description='No')
    buttons = HBox([yes,no])
    display(buttons)


def notebook_ui():
    print("cif or pdb?\n")
    pdb = widgets.Button(description='pdb')
    cif = widgets.Button(description='cifs')
    buttons = HBox([pdb,cif])
    display(buttons)
    pdb.on_click(pdbbutton)
    cif.on_click(cifbutton)
  