### Command-line python script to compute TC/L interaction probability of compound
### Version 1.0 (27.11.22): https://github.com/juppifluppi/tclint

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from scopy.ScoPretreat import pretreat
import scopy.ScoDruglikeness
from dimorphite_dl import DimorphiteDL
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# load training and test set data

comp = 1

try:

    df  = pd.read_csv("trainvalues.csv")
    x = df['rd_logD']
    y = df['rd_MR']

    df2  = pd.read_csv("testvalues.csv")
    x2 = df2['rd_logD']
    y2 = df2['rd_MR']

    file1 = open('smiles_train.csv', 'r')
    fil = file1.readlines()

    file2 = open('smiles_test.csv', 'r')
    fil2 = file2.readlines()

except:
    print("Training/test data sets not found, continuing without comparison.")
    comp = 0
    sys.exit()

# convert to rdkit mols    

if comp == 1:

    o=[]
    for lines in fil:
        o.append(Chem.MolFromSmiles(lines))
    
    ox=[]
    for lines in fil2:
        ox.append(Chem.MolFromSmiles(lines))

# protonate and pretreat given SMILES, compute descriptors via rdkit/scopy, calculate probability    
    
try:

    SMI = str(sys.argv[1])
    
    dimorphite_dl = DimorphiteDL(
        min_ph = 6.4,
        max_ph = 6.6,
        max_variants = 1,
        label_states = False,
        pka_precision = 0.1
    )
    SMI = str(dimorphite_dl.protonate(SMI)[0])
    
    mol = Chem.MolFromSmiles(SMI)
    sdm = pretreat.StandardizeMol()
    mol = sdm.disconnect_metals(mol)
     
    logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(mol)
    mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(mol)
    
    tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
    tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333
    
    tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )
    
    print("TC/L interaction probability model")
    print("---------------------------------")
    print("Compound properties:")  
    print("logD: " + str(round(logd,2)))
    print("CrippenMR: " + str(round(mr,2)))
    print("TC/L interaction probability: " + str(round(tcl3*100,2)) + " %")
    print("---------------------------------")  

except:
    print("Something is wrong with your SMILES code.")
    sys.exit()

# copmute ecfp_4 fingerprints to calculate SDC metrics    

if comp == 1:

    fp1 = AllChem.GetMorganFingerprint(mol, 2)

    g=[]

    for molx in o:
        fp2 = AllChem.GetMorganFingerprint(molx, 2)
        Tan = DataStructs.TanimotoSimilarity(fp1,fp2)
        try:
            scd = 2.718281828459045 ** ((-3 * Tan)/(1 - Tan))
            g.append(scd)
        except:
            pass

    li=[]
    
    for k in ox:
        lkx=[]
        fp1 = AllChem.GetMorganFingerprint(k, 2)
        for molx in o:
            fp2 = AllChem.GetMorganFingerprint(molx, 2)
            Tan = DataStructs.TanimotoSimilarity(fp1,fp2)
            try:
                scd = 2.718281828459045 ** ((-3 * Tan)/(1 - Tan))
                lkx.append(scd)
            except:
                pass
        li.append(np.sum(lkx))
 
    lit=[] 

    for k in o:
        lkx=[]
        fp1 = AllChem.GetMorganFingerprint(k, 2)
        for molx in o:
            fp2 = AllChem.GetMorganFingerprint(molx, 2)
            Tan = DataStructs.TanimotoSimilarity(fp1,fp2)
            try:
                scd = 2.718281828459045 ** ((-3 * Tan)/(1 - Tan))
                lkx.append(scd)
            except:
                pass
        lit.append(np.sum(lkx)) 

    # output values    
    
    print("Compound vs modeling sets:")
    print("Training set logD: "+str(round(np.min(x),2))+" - "+str(round(np.max(x),2))+" (Mean: "+str(round(np.mean(x),2))+"; SD: "+str(round(np.std(x),2))+")")       
    print("Training set CrippenMR: "+str(round(np.min(y),2))+" - "+str(round(np.max(y),2))+" (Mean: "+str(round(np.mean(y),2))+"; SD: "+str(round(np.std(y),2))+")")  
    print("Validation set logD: "+str(round(np.min(x2),2))+" - "+str(round(np.max(x2),2))+" (Mean: "+str(round(np.mean(x2),2))+"; SD: "+str(round(np.std(x2),2))+")")       
    print("Validation set CrippenMR: "+str(round(np.min(y2),2))+" - "+str(round(np.max(y2),2))+" (Mean: "+str(round(np.mean(y2),2))+"; SD: "+str(round(np.std(y2),2))+")")  

    print("---------------------------------")         
    print("SDC applicability domain metrics:")
    print("Training set SDC: "+str(round(np.min(lit),2))+" - "+str(round(np.max(lit),2))+" (Mean: "+str(round(np.mean(lit),2))+"; SD: "+str(round(np.std(lit),2))+")")       
    print("Validation set SDC: "+str(round(np.min(li),2))+" - "+str(round(np.max(li),2))+" (Mean: "+str(round(np.mean(li),2))+"; SD: "+str(round(np.std(li),2))+")")     
    print("Compound SDC: "+str(round(np.sum(g),2)))
    print("---------------------------------")   
print("Version: 1.0 (27.11.22)")  
print("Web version: https://tclint.streamlit.app/")
print("More info can be found on: https://github.com/juppifluppi/tclint")   
