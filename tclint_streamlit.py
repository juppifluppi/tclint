### Streamlit web app to compute TC/L interaction probability of compound
### Josef Kehrein
### Version 1.4 (29.12.22): https://github.com/juppifluppi/tclint

# load modules and define protonation settings

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from scopy.ScoPretreat import pretreat
import scopy.ScoDruglikeness
from dimorphite_dl import DimorphiteDL
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

dimorphite_dl = DimorphiteDL(
    min_ph = 6.4,
    max_ph = 6.6,
    max_variants = 1,
    label_states = False,
    pka_precision = 0.1
)

# description of model

st.title('TC/L interaction probability model')
st.caption("""Input a [molecule SMILES code](https://pubchem.ncbi.nlm.nih.gov/edit3/index.html).
A taurocholate/lecithin interaction probability is computed at pH 6.5, based on two descriptors: logD and CrippenMR.
The model is based on [our previous study](https://doi.org/10.1021/acs.molpharmaceut.2c00227)
but was built with descriptors from rdkit/scopy instead of MOE/PaDEL, using logD for pH 7.4 instead of 7.0.
For the same traning/validation set similar statistics are retrieved (balanced accuracy: 0.86/0.83, AUC: 0.93/0.93).""")

# load values for training set used for model building (SMILES, logD, MR, class)

train_SMI = ['Clc1c(F)ccc(Nc2ncnc3c2cc(NC(=O)C=CC[NH+](C)C)c(OC2COCC2)c3)c1', 'O=C1NC=Nc2[nH]ncc12', 'FC(F)(F)C(C)(C)c1nccc(-c2c(C)nc(NC(=O)N3C(C(=O)N)CCC3)s2)c1', 'Brc1c(N)c(C[NH2+]C2CCC(O)CC2)cc(Br)c1', 'Ic1c(OCC[NH+](CC)CC)c(I)cc(C(=O)c2c(CCCC)oc3c2cccc3)c1', '[NH+](CCC=C1c2c(cccc2)CCc2c1cccc2)(C)C', 'Clc1c(C2C(C(=O)OCC)=C(COCC[NH3+])NC(C)=C2C(=O)OC)cccc1', 'O=C(NC1C(=O)N2C(C(=O)[O-])C(C)(C)SC12)C([NH3+])c1ccc(O)cc1', 'ClC(F)(F)Oc1ccc(NC(=O)c2cc(c(N3CC(O)CC3)nc2)-c2[nH]ncc2)cc1', 'O=C(Oc1c(C(=O)[O-])cccc1)C', 'O=C(OC)NC(C(=O)NN(CC(O)C(NC(=O)C(NC(=O)OC)C(C)(C)C)Cc1ccccc1)Cc1ccc(-c2ncccc2)cc1)C(C)(C)C', 'O=C(N)Cc1ccc(OCC(O)C[NH2+]C(C)C)cc1', 'O=C(OC1CC2[NH+](C)C(C1)CC2)C(CO)c1ccccc1', 'O=C(OCC)c1ccc(N)cc1', 'O=C(Oc1ccc(C(c2ncccc2)c2ccc(OC(=O)C)cc2)cc1)C', 'O=C(NC)c1nccc(Oc2cc3sc(NC4C(O)CCCC4)nc3cc2)c1', 'Brc1c(N)c(C[NH+](C)C2CCCCC2)cc(Br)c1', 'O=C1N(CCCC[NH+]2CCN(c3ncccn3)CC2)C(=O)CC2(C1)CCCC2', 'O=C([O-])c1c2n(Cc3ccc(-c4c(-c5nnn[n-]5)cccc4)cc3)c(OCC)nc2ccc1', 'Fc1c(C(=O)NC)ccc(C2=Nn3c(Cc4cc5c(nccc5)cc4)cnc3N=C2)c1', 'O=C(N)N1c2c(cccc2)C=Cc2c1cccc2', 'S=C1N(C(=O)OCC)C=CN1C', 'S(=O)(=O)(N)c1ccc(-n2c(-c3ccc(C)cc3)cc(C(F)(F)F)n2)cc1', 'Clc1c(Nc2c(S(=O)(=O)C(C)C)cccc2)nc(Nc2c(OC(C)C)cc(c(C)c2)C2CC[NH2+]CC2)nc1', 'Clc1ccc(C(N2CC[NH+](CCOCC(=O)[O-])CC2)c2ccccc2)cc1', 'O(C)c1cc2c(C(O)C3[NH+]4CC(C=C)C(C3)CC4)ccnc2cc1', 'ClC(Cl)C(=O)NC(C(O)c1ccc([N+](=O)[O-])cc1)CO', 'Clc1cc2[nH+]ccc(NC(CCC[NH+](CC)CC)C)c2cc1', 'Clc1cc2N(CCC[NH+](C)C)c3c(Sc2cc1)cccc3', 'Clc1cc2C(=CCC[NH+](C)C)c3c(Sc2cc1)cccc3', 'C(=Cc1ccccc1)C[NH+]1CCN(C(c2ccccc2)c2ccccc2)CC1', 'Clc1c([NH+]=C2NCCN2)c(Cl)ccc1', 'Clc1c(C(C(=O)OC)N2Cc3c(scc3)CC2)cccc1', 'O=C1N(C)C(=O)c2n(C)cnc2N1C', 'O=C(CO)C1(O)C2(C)C(C3C(C(=O)C2)C2(C)C(=CC(=O)CC2)CC3)CC1', 'OC1(C#C)C2(C)C(C3C(C4(C)C(=Cc5oncc5C4)CC3)CC2)CC1', 'FC12C3(C)C(=CC(=O)C=C3)CCC1C1C(C)(C(O)(C(=O)CO)C(C)C1)CC2O', 'O(C)c1cc2c(cc1)CC1[NH+](C)CCC32C1CCCC3', 'Clc1c(Nc2c(CC(=O)[O-])cccc2)c(Cl)ccc1', 'O=C(OC1C(=O)N(CC[NH+](C)C)c2c(SC1c1ccc(OC)cc1)cccc2)C', 'Clc1cc2NC(=O)N(C3CC[NH+](CCCN4C(=O)Nc5c4cccc5)CC3)c2cc1', 'O=C([O-])c1cc(-c2c(O)c(N=NC3C(=O)N(c4cc(C)c(C)cc4)N=C3C)ccc2)ccc1', 'Clc1ccc(C(=O)c2ccc(OC(C(=O)OC(C)C)(C)C)cc2)cc1', 'O=C([O-])C(NC(=O)c1ccc(NCc2nc3C(=O)N=C(N)Nc3nc2)cc1)CCC(=O)[O-]', 'Clc1c(S(=O)(=O)N)cc(C(=O)[O-])c(NCc2occc2)c1', 'Clc1cc(C(=O)NCCc2ccc(S(=O)([O-])=NC(=O)NC3CCCCC3)cc2)c(OC)cc1', 'S(=O)(=O)(NC(=O)NC1CCC(C)CC1)c1ccc(CCC(=O)NN2C(=O)C(CC)=C(C)C2)cc1', 'Clc1ccc(C2(O)CC[NH+](CCCC(=O)c3ccc(F)cc3)CC2)cc1', 'Clc1c(S(=O)(=O)N)cc2S(=O)(=O)NCNc2c1', 'O=C(CO)C1(O)C2(C)C(C3C(C4(C)C(=CC(=O)CC4)CC3)C(O)C2)CC1', 'O=C(Nc1cc(Nc2nc(-c3cnccc3)ccn2)c(C)cc1)c1ccc(CN2CC[NH+](C)CC2)cc1', '[NH+](CCCN1c2c(cccc2)CCc2c1cccc2)(C)C', 'Clc1ccc(C(=O)n2c(C)c(CC(=O)[O-])c3c2ccc(OC)c3)cc1', 'O=C(NN)c1ccncc1', 'OC(C[NH2+]C(C)C)c1cc(O)c(O)cc1', 'Clc1c(C2(Cn3ncnc3)OC(COc3ccc(N4CCN(c5ccc(N6C(=O)N(C(CC)C)N=C6)cc5)CC4)cc3)CO2)ccc(Cl)c1', 'O=C1OC2CC3(OC(C(CC)C)C(C)CC3)OC(CC=C(C)C(OC3OC(C)C(OC4OC(C)C(O)C(OC)C4)C(OC)C3)C(C)C=CC=C3C4(O)C(C(O)C(C)=CC14)OC3)C2', 'O=C([O-])C(C)c1cc(C(=O)c2ccccc2)ccc1', 'Clc1c(OCc2cc(F)ccc2)ccc(Nc2ncnc3c2cc(-c2oc(C[NH2+]CCS(=O)(=O)C)cc2)cc3)c1', 'O=C(Nc1c(C)cccc1C)C[NH+](CC)CC', 'Clc1cc2c(C(=C3CCN(C(=O)OCC)CC3)c3ncccc3CC2)cc1', 'FC(F)(C)c1c(-c2c(Oc3ccc(C=CC(=O)[O-])cc3)c3c(s2)cc(O)cc3)ccc(F)c1', 'FC(F)(F)c1nccc(C(=O)Nc2cc(c(C)cc2)-c2cc(OCCO)nc(N3CCOCC3)c2)c1', '[NH2+](CCCC12c3c(cccc3)C(c3c1cccc3)CC2)C', 'FC(F)(F)c1nc2c(C(F)(F)F)cccc2c(C(O)C2[NH2+]CCCC2)c1', 'O=C1C(C)=CC(=O)c2c1cccc2', 'O=C(Nc1c(C)cccc1C)C1[NH+](C)CCCC1', 'S(=O)(=O)([O-])CN(C)C=1C(=O)N(N(C)C=1C)c1ccccc1', '[NH2+]=C(N=C(N)N)N(C)C', 'Clc1c(N)cc(OC)c(C(=O)NCC[NH+](CC)CC)c1', 'O(CCc1ccc(OCC(O)C[NH2+]C(C)C)cc1)C', 'O=[N+]([O-])c1n(CCO)c(C)nc1', 'O=C(N)C=1C(=O)C2(O)C(O)=C3C(=O)c4c(O)ccc(N(C)C)c4CC3CC2C([NH+](C)C)C=1[O-]', 'Clc1c([N-]C(=O)c2c(O)ccc(Cl)c2)ccc([N+](=O)[O-])c1', 'O=[N+]([O-])c1c(C2C(C(=O)OC)=C(C)NC(C)=C2C(=O)OC)cccc1', 'O=[N+]([O-])c1oc(C=NN2C(=O)NC(=O)C2)cc1', 'O=C([O-])c1c(O)cc(N)cc1', 'O(C)c1c(OC)ccc(Cc2[nH+]ccc3c2cc(OC)c(OC)c3)c1', 'O=C(Nc1ccc(O)cc1)C', 'SC(C([NH3+])C(=O)[O-])(C)C', 'O=C(NC1C(=O)N2C(C(=O)[O-])C(C)(C)SC12)COc1ccccc1', 'Clc1cc2N(CCC[NH+]3CCN(CCO)CC3)c3c(Sc2cc1)cccc3', 'O=C1N(N(C)C(C)=C1)c1ccccc1', 'O=C1C(C(CC)c2ccccc2)=C([O-])c2c(O1)cccc2', 'O=C1C(c2ccccc2)(c2ccccc2)NC(=O)N1', 'Fc1c(C2(Cn3ncnc3)OCC(COc3ccc(N4CCN(c5ccc(N6C(=O)N(C(C(O)C)CC)N=C6)cc5)CC4)cc3)C2)ccc(F)c1', 'O=C(OCC(=O)C1(O)C2(C)C(C3C(C4(C)C(=CC(=O)C=C4)CC3)C(O)C2)CC1)C', 'O=C(OCC[NH+](CC)CC)c1ccc(N)cc1', '[NH+](C(CN1c2c(Sc3c1cccc3)cccc2)C)(C)C', 'O(CC(O)C[NH2+]C(C)C)c1c2c(ccc1)cccc2', 'O=C1N(N(C)C(C)=C1CCC)c1ccccc1', 'O=C1N(C)C(=O)c2n(CC(O)C)cnc2N1C', 'O(CCO)CCN1CCN(C2=[NH+]c3c(Sc4c2cccc4)cccc3)CC1', 'O=C(OCC)C(NC(C(=O)N1C(C(=O)[O-])CC2C1CCC2)C)CCc1ccccc1', 'S(CCNC(NC)=C[N+](=O)[O-])Cc1oc(C[NH+](C)C)cc1', 'O=C([O-])c1c(O)cccc1', 'S(C(=O)C)C1C2C3C(C)(C4(OC(=O)CC4)CC3)CCC2C2(C)C(=CC(=O)CC2)C1', 'S(=O)([O-])(=Nc1noc(C)c1)c1ccc(N)cc1', 'S(=O)(=O)(NC)Cc1cc2c(CC[NH+](C)C)c[nH]c2cc1', 'Brc1c(N)nc(-n2nccc2)nc1-n1nccc1', 'O=C(OCC[NH+](C)C)c1ccc(NCCCC)cc1', 'O=C(N)C=1C(=O)C2(O)C(O)=C3C(=O)c4c(O)cccc4C(O)(C)C3CC2C([NH+](C)C)C=1[O-]', 'O=C1N(C)C(=O)c2nc[nH]c2N1C', 'S(=O)([O-])(=NC(=O)NCCCC)c1ccc(C)cc1', 'Ic1cc(F)c(NC=2N(C)C(=O)C(C)=C3N(c4cc(NC(=O)C)ccc4)C(=O)N(C(=O)C=23)C2CC2)cc1', 'O(C)c1c(OC)cc(Cc2c(N)nc(N)[nH+]c2)cc1OC', 'O=C([O-])C(CCC)CCC', 'O(C)c1c(OC)ccc(C(C#N)(C(C)C)CCC[NH+](CCc2cc(OC)c(OC)cc2)C)c1']
train_logD = [2.72, 0.15, 2.21, 2.42, 5.34, 2.15, 0.69, -2.96, 2.69, -0.82, 3.53, 0.45, 1.03, 1.49, 3.06, 2.17, 3.17, 0.46, 1.26, 3.2, 2.95, 1.82, 2.6, 3.83, -0.62, 2.07, 1.11, 2.34, 3.14, 2.87, 2.51, 1.31, 2.04, 0.02, 1.75, 4.57, 1.82, 2.0, 1.48, 1.81, 1.63, 0.41, 3.75, -3.43, 0.05, 1.93, 1.38, 2.63, 0.11, 1.76, 1.78, 2.41, 1.13, 0.52, 0.32, 4.62, 5.61, 0.53, 4.19, 1.25, 2.78, 3.39, 3.22, 2.55, 2.69, 1.78, 1.29, -0.64, -3.46, 1.35, 1.2, 0.82, -1.55, 3.18, 1.3, 0.79, -1.06, 3.17, 1.2, -1.74, -1.48, 2.03, 1.41, 3.45, 1.58, 4.01, 1.76, 1.03, 2.71, 1.8, 2.07, -0.16, 1.03, -1.24, -0.17, -0.88, 3.83, 0.68, 0.39, 1.79, 1.35, -2.1, -0.19, 0.52, 2.71, 1.67, 0.06, 2.82]
train_MR = [129.53, 34.51, 106.23, 80.15, 142.54, 89.67, 103.07, 85.89, 110.48, 42.08, 194.1, 72.77, 79.08, 46.81, 100.68, 107.17, 83.72, 105.02, 119.38, 114.12, 73.53, 46.66, 90.11, 152.95, 101.82, 94.15, 72.57, 93.95, 91.32, 91.39, 116.46, 58.1, 83.77, 51.2, 94.14, 96.35, 99.94, 80.49, 74.9, 113.62, 118.58, 122.69, 97.26, 106.64, 73.19, 127.35, 127.12, 100.11, 61.64, 95.14, 145.14, 89.59, 93.12, 36.11, 56.66, 188.18, 227.4, 69.74, 152.91, 71.46, 106.07, 119.01, 127.22, 86.69, 81.14, 48.86, 73.94, 76.86, 33.83, 81.66, 75.45, 40.7, 111.87, 79.18, 87.19, 53.42, 36.85, 97.2, 42.41, 34.75, 84.09, 111.85, 55.75, 82.72, 70.56, 187.31, 104.6, 68.05, 86.29, 77.38, 69.74, 61.9, 106.78, 107.08, 81.1, 32.44, 112.21, 63.22, 81.36, 68.39, 77.8, 103.41, 46.58, 70.95, 148.76, 79.76, 38.31, 130.78]
train_class = [1,0,1,1,1,1,1,0,1,0,1,0,0,1,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,1,1,0,1,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1]

# load values for validation set used for model evaluation (SMILES, logD, MR, class)

test_SMI = ['S(=O)(=O)(N)c1sc(NC(=O)C)nn1', 'SCC(C(=O)N1C(C(=O)[O-])CCC1)C', 'O(C)c1cc2c(C(O)C3[NH+]4CC(C=C)C(C3)CC4)ccnc2cc1', '[NH+](CCC1=C(C(C)c2ncccc2)c2c(cccc2)C1)(C)C', 'O(C(c1ccccc1)c1ccccc1)CC[NH+](C)C', '[NH+](CCC=C1c2c(OCc3c1cccc3)cccc2)(C)C', 'O=C(N)C=1C(=O)C2(O)C(O)=C3C(=O)c4c(O)cccc4C(C)C3C(O)C2C([NH+](C)C)C=1[O-]', 'O=C1N(C)C(=O)c2n(CCO)cnc2N1C', 'O=C([O-])C(C)c1ccc(CC(C)C)cc1', 'Clc1ccc(C2(O)CC[NH+](CCC(C(=O)N(C)C)(c3ccccc3)c3ccccc3)CC2)cc1', 'O=C(Oc1c(C)c(C)c(OCC(O)C[NH2+]C(C)C)cc1C)C', 'FC(F)(F)c1cc(NC(=O)c2cc(Nc3nc(-c4cnccc4)ccn3)c(C)cc2)cc(-n2cc(C)nc2)c1', '[SH0](=O)(Cc1c(C)c(OC)c(C)cn1)c1[nH]c2c(n1)cc(OC)cc2', 'O=C(Nc1ccc(OCC)cc1)C', 'O=C(Nc1ncccc1)C=1N(C)S(=O)(=O)c2c(C=1[O-])cccc2', 'OCc1c(O)c(C)ncc1CO', 'OC(C[NH2+]C(C)(C)C)c1cc(CO)c(O)cc1', 'Nc1c(-c2ccccc2)nc2c(N)nc(N)nc2n1']
test_logD = [-0.14, -1.53, 2.07, 2.77, 1.96, 1.91, -2.01, -0.46, 0.56, 1.97, 1.1, 4.34, 1.5, 1.79, 0.08, 0.33, 0.39, 1.57]
test_MR = [45.59, 51.88, 94.15, 92.01, 78.36, 86.59, 103.64, 57.3, 58.41, 136.25, 85.39, 140.98, 93.02, 51.91, 82.93, 42.48, 65.53, 73.8]
test_class = [0,0,1,1,1,1,1,0,0,1,0,1,1,1,0,0,0,0]

# load values for threshold line (equalling 50% probability)

thresh_x = [-4.0,-3.8,-3.6,-3.4,-3.2,-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8]
thresh_y = [155,155,150,150,145,145,140,135,135,130,130,125,125,120,120,115,110,110,105,105,100,100,95,90,90,85,85,80,80,75,70,70,65,65,60,60,55,50,50,45,45,40,40,35,35,30,25,25,20,20]

# protonate training set SMILES via dimorphite-dl at pH 6.5 and convert to rdkit mol, then calculate probability

train_mols = []
train_prob = []
for lines in train_SMI:
    SMI = str(dimorphite_dl.protonate(lines)[0])
    mol = Chem.MolFromSmiles(SMI)
    sdm = pretreat.StandardizeMol()
    mol = sdm.disconnect_metals(mol)
    train_mols.append(mol)    
    logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(mol)
    mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(mol)
    tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
    tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333 
    tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )
    train_prob.append(tcl3*100)

# do the same for the validation set

test_mols = []
test_prob = []
for lines in test_SMI:
    SMI = str(dimorphite_dl.protonate(lines)[0])
    mol = Chem.MolFromSmiles(SMI)
    sdm = pretreat.StandardizeMol()
    mol = sdm.disconnect_metals(mol)
    test_mols.append(mol)    
    logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(mol)
    mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(mol)
    tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
    tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333 
    tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )
    test_prob.append(tcl3*100)

# do the same for the input SMILES of the user, output molecular descriptors and probability

try:
    SMI = st.text_input('Enter SMILES code', 'CC(C)NCC(COC1=CC=C(C=C1)CCOC)O')  
    SMI = str(dimorphite_dl.protonate(SMI)[0])    
    mol = Chem.MolFromSmiles(SMI)
    sdm = pretreat.StandardizeMol()
    mol = sdm.disconnect_metals(mol)    
    m = Chem.MolFromSmiles(SMI)
    im = Draw.MolToImage(m,fitImage=True)    
    logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(mol)
    mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(mol)    
    tcl1 = ( ( logd - 1.510648) / 1.708574 ) * 1.706694
    tcl2 = ( ( mr - 90.62889 ) / 35.36033 ) * 2.4925333    
    tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.9872289 + tcl1 + tcl2 ) ) ) )    
    st.image(im)
    st.write("logD: " + str(round(logd,2)))
    st.write("CrippenMR: " + str(round(mr,2)))
    st.write("TC/L interaction probability: " + str(int(round(tcl3*100,2))) + " %")
    
# output error if SMILES cannot be recognized

except:
    st.write("Something is wrong with your SMILES code.")
    st.stop()
    
# description of additional plot and SDC metric

st.caption("""The following plot shows the properties in relation to the modeling sets. A [SDC metric](https://doi.org/10.1021/acs.jcim.8b00597) (sum of tanimoto distance-weighted contributions) evaluates structural similarity to the training set. Higher SDC values and / or large distances to the training set in the plot can indicate less reliable predictions.""")

# compute tanimoto ecfp_4 fingerprints (approximately radius = 2 of MorganFingerprint) to calculate SDC metrics

# compare user compound to training set

fp1 = AllChem.GetMorganFingerprint(mol, 2)
compound_sdc=[]
for k in train_mols:
    fp2 = AllChem.GetMorganFingerprint(k, 2)
    Tan = DataStructs.TanimotoSimilarity(fp1,fp2)
    try:
        sdc = 2.718281828459045 ** ((-3 * Tan)/(1 - Tan))
        compound_sdc.append(sdc)
    except:
        pass
 
# compare training set compounds between each other 
 
train_sdc=[]
for k in train_mols:
    values=[]
    fp1 = AllChem.GetMorganFingerprint(k, 2)
    for molx in train_mols:
        fp2 = AllChem.GetMorganFingerprint(molx, 2)
        Tan = DataStructs.TanimotoSimilarity(fp1,fp2)
        try:
            sdc = 2.718281828459045 ** ((-3 * Tan)/(1 - Tan))
            values.append(sdc)
        except:
            pass
    train_sdc.append(np.sum(values)) 

# plot values of training and validation sets
    
fig=plt.figure()
ax=fig.add_axes([0,0,1,1])
l=ax.scatter(train_logD, train_MR, color='b',alpha=0.5,s=50)
p=ax.scatter(test_logD, test_MR, color='r',alpha=0.5,s=50)
ax.set_xlabel('logD')
ax.set_ylabel('CrippenMR')
ax.set_title('Compound vs. modeling sets')

# add additional markers to training set for missclassified instances and > 50 % probability

for j in range(0,len(train_prob)):
    if train_prob[j] >= 50:
        m=ax.scatter(train_logD[j],train_MR[j], color='b',alpha=0.1,s=160)
        if train_class[j] == 0:
            n=ax.scatter(train_logD[j],train_MR[j], color='y',alpha=1,s=10,zorder=2)
    if train_prob[j] < 50:
        if train_class[j] == 1:
            h=ax.scatter(train_logD[j],train_MR[j], color='y',alpha=1,s=10,zorder=2)
            
# do the same for the validation set
            
for j in range(0,len(test_prob)):
    if test_prob[j] >= 50:
        w=ax.scatter(test_logD[j],test_MR[j], color='r',alpha=0.1,s=160)
        if test_class[j] == 0:
            r=ax.scatter(test_logD[j],test_MR[j], color='y',alpha=1,s=10,zorder=2)
    if test_prob[j] < 50:
        if test_class[j] == 1:
            s=ax.scatter(test_logD[j],test_MR[j], color='y',alpha=1,s=10,zorder=2)

# plot user compound
            
o=ax.scatter(logd, mr, color='lime',alpha=1,s=50,marker="D",zorder=2)            
            
# use calculated values for threshold line to fit regression line

b, a = np.polyfit(thresh_x, thresh_y, deg=1)
xseq = np.linspace(-4, 5.5, num=2)
ax.plot(xseq, a + b * xseq, color="grey", linewidth=10, alpha=0.1,zorder=1)
knn=ax.scatter(xseq, a + b * xseq, color="grey", s=30, alpha=0.1,edgecolors="grey", marker="s",zorder=1)

# add legend to plot

ax.legend((l,m,p,w,o,n,knn),("Training set", "≥50% probability", "Validation set", "≥50% probability", "Compound", "Misclassified", "50% threshold"), ncol=1)
plt.show()
st.pyplot(fig)

# output SDC metrics for training set and compound

st.write("SDC applicability domain metrics:")
st.write("Training set: "+str(round(np.min(train_sdc),2))+" - "+str(round(np.max(train_sdc),2))+" (Mean: "+str(round(np.mean(train_sdc),2))+"; SD: "+str(round(np.std(train_sdc),2))+")")       
st.write("Compound: "+str(round(np.sum(compound_sdc),2)))

# reference

st.caption("Version 1.4 (29.12.22). Visit [github](https://github.com/juppifluppi/tclint) for more information and a downloadable version.")
