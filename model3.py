from rdkit import Chem
from rdkit.Chem import Draw
from padelpy import from_smiles
from scopy.ScoPretreat import pretreat
from dimorphite_dl import DimorphiteDL
import streamlit as st

st.header('Lipid interaction probability model')
#st.caption("""Input a SMILES code of your molecule of choice (use e.g. https://pubchem.ncbi.nlm.nih.gov/edit3/index.html).
#A probability for interaction with taurocholate/lecithin is computed for the compound at pH 6.5, based on two descriptors: logD and CrippenMR.
#The model is based on Mol. Pharmaceutics 2022, 19, 2868−2876 (https://doi.org/10.1021/acs.molpharmaceut.2c00227),
#but descriptors are computed via rdkit/scopy instead of MOE/PaDEL, and logD for pH 7.4 instead of 7.0 is used.""")



SMI = st.text_input('Enter SMILES of drug molecule', 'CC(C)NCC(COC1=CC=C(C=C1)CCOC)O')


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
    
descriptors = from_smiles(SMI)
    
ve2 = float(descriptors["VE2_Dt"])
sw = float(descriptors["SwHBa"])
lgp = float(descriptors["CrippenLogP"])
lgp = ( lgp * 0.8469470 ) + 0.2950044
hb = float(descriptors["ETA_Beta_ns"])
hb = ( hb * 0.7306608 ) - 0.1405840
    
lipid1 = ( ( ve2 - 0.01 ) / 0.01 ) * 0.782 * (-1)
lipid2 = ( ( sw - 12.85 ) / 8.96 ) * 0.624   
lipid3 = ( ( lgp - 2.10 ) / 1.99 ) * 0.474 
lipid4 = ( ( hb - 7.12 ) / 3.35 ) * 0.036  
   
lipid5 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( lipid1 + lipid2 + lipid3 + lipid4 - 8.088 ) ) ) )
    
st.image(im)
st.text("VE2_Dt: " + str(round(lipid1,2)))
st.text("SwHBa: " + str(round(lipid2,2)))
st.text("SlogP: " + str(round(lipid3,2)))
st.text("h_log_pbo: " + str(round(lipid4,2)))
st.text("Lipid interaction probability: " + str(round(lipid5,2)))
