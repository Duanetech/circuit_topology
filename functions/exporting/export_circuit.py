import pandas as pd

def export_circuit(circlist):
    df = pd.DataFrame(circlist,columns = ['protid','segnums','meanlength','segends'])
    df.to_csv('results/circuit/circuitlist.csv',index=False)
    print('Succesfully exported circuitlist.csv')