from rdkit import Chem
from rdkit.Chem import Draw
from scopy.ScoPretreat import pretreat
import scopy.ScoDruglikeness
from dimorphite_dl import DimorphiteDL
import streamlit as st

import pickle
import datamaker
import streamlit as st
st.title('Streamlit + RDKit :rocket:')
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from io import BytesIO
from functools import partial
from PIL import Image
from rdkit.Chem.Draw import rdDepictor
rdDepictor.SetPreferCoordGen(True)

st.header('TC/L interaction probability model')
st.caption("""Input a SMILES code of your molecule of choice.
A probability for interaction with taurocholate/lecithin is computed based on two descriptors, logD and CrippenMR.
The model is based on results from Mol. Pharmaceutics 2022, 19, 2868−2876, but descriptors are computed via rdkit/scopy instead of MOE/PaDEL.""")

import streamlit.components.v1 as components

_RELEASE = False

if not _RELEASE:
    _component_func = components.declare_component(
            "chemstreamlit",
            url="http://<public_ip>:3001"
            )
else:
    parent_dir = os.path.dirname(os.path.abspath(__file__))
    build_dir = os.path.join(parent_dir, "frontend/build")
    _component_func = components.declare_component("my_component", path=build_dir)

def my_component():
    component_value = _component_func()
    return component_value

SMI = my_component()

#SMI = st.text_input('Enter Canonical SMILES of drug', 'CC(C)NCC(COC1=CC=C(C=C1)CCOC)O')

dimorphite_dl = DimorphiteDL(
    min_ph=6.5,
    max_ph=6.5,
    max_variants=1,
    label_states=False,
    pka_precision=0.5
)
SMI = str(dimorphite_dl.protonate(SMI)[0])

mol = Chem.MolFromSmiles(SMI)
sdm = pretreat.StandardizeMol()
mol = sdm.disconnect_metals(mol)

m = Chem.MolFromSmiles(SMI)
im=Draw.MolToImage(m)

logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(mol)
mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(mol)

tcl1 = ( ( logd - 2.009365 ) / 1.878026 ) * 1.1626256
tcl2 = ( ( mr - 90.80861 ) / 35.30108 ) * 1.9764294

tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.7999132 + tcl1 + tcl2 ) ) ) )

st.image(im)
st.text("logD: " + str(round(logd,2)))
st.text("CrippenMR: " + str(round(mr,2)))
st.text("TC/L interaction probability: " + str(round(tcl3,2)))
