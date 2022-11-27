from rdkit import Chem
from rdkit.Chem import Draw
from scopy.ScoPretreat import pretreat
import scopy.ScoDruglikeness
from dimorphite_dl import DimorphiteDL
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem

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

o=[]
for lines in fil:
    o.append(Chem.MolFromSmiles(lines))
    
ox=[]
for lines in fil2:
    ox.append(Chem.MolFromSmiles(lines))

st.title('TC/L interaction probability model')
st.caption("""Input a SMILES code of a molecule (use e.g. https://pubchem.ncbi.nlm.nih.gov/edit3/index.html).
A probability for interaction with taurocholate/lecithin is computed at pH 6.5, based on two descriptors: logD and CrippenMR.
The model is inspired by Mol. Pharmaceutics 2022, 19, 2868−2876 (https://doi.org/10.1021/acs.molpharmaceut.2c00227)
but was built with descriptors from rdkit/scopy instead of MOE/PaDEL, using logD for pH 7.4 instead of 7.0.
For the same traning/validation set similar statistics are retrieved (balanced accuracy: 0.86/0.83, AUC: 0.93/0.93).""")


try:

    SMI = st.text_input('Enter SMILES', 'CC(C)NCC(COC1=CC=C(C=C1)CCOC)O')
    
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
    
    m = Chem.MolFromSmiles(SMI)
    im = Draw.MolToImage(m)
    
    logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(mol)
    mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(mol)
    
    tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
    tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333
    
    tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )
    
    st.image(im)
    st.write("logD: " + str(round(logd,2)))
    st.write("CrippenMR: " + str(round(mr,2)))
    st.write("TC/L interaction probability: " + str(round(tcl3,2)))

except:
    pass

st.caption("""A scatter plot shows the properties of the compound in relation to the training and the validation set. A SDC applicability
domain metric based on summed tanimoto similarities to the training set molecules is given underneath (J. Chem. Inf. Model. 2019, 59, 
181−189). Higher SDC values and / or large distances to the training set in the plot can indicate less reliable predictions.""")

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
         
fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
ax.scatter(x, y, color='b',alpha=0.5)
ax.scatter(x2, y2, color='r',alpha=0.5)
ax.scatter(logd, mr, color='g',alpha=1,s=150)
ax.set_xlabel('logD')
ax.set_ylabel('CrippenMR')
ax.set_title('Compound vs. modeling sets')

l=ax.scatter(x, y, color='b',alpha=0.5)
p=ax.scatter(x2, y2, color='r',alpha=0.5)
o=ax.scatter(logd, mr, color='g',alpha=1)
ax.legend((l,p,o),("Training set", "Validation set", "Compound"))
plt.show()

st.pyplot(fig)

st.write("SDC applicability domain metrics:")
st.write("Training set: "+str(round(np.min(lit),2))+" - "+str(round(np.max(lit),2))+" (Mean: "+str(round(np.mean(lit),2))+"; SD: "+str(round(np.std(lit),2))+")")       
st.write("Validation set: "+str(round(np.min(li),2))+" - "+str(round(np.max(li),2))+" (Mean: "+str(round(np.mean(li),2))+"; SD: "+str(round(np.std(li),2))+")")     
st.write("Compound: "+str(round(np.sum(g),2)))


