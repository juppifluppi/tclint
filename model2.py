from rdkit import Chem
from rdkit.Chem import Draw
from scopy.ScoPretreat import pretreat
import scopy.ScoDruglikeness
from dimorphite_dl import DimorphiteDL
import streamlit as st

SMI='CC(C)(C)C1=CC(=CC(=C1O)C(C)(C)C)SC(C)(C)SC2=CC(=C(C(=C2)C(C)(C)C)O)C(C)(C)C'

dimorphite_dl = DimorphiteDL(
    min_ph=6.5,
    max_ph=6.5,
    max_variants=1,
    label_states=False,
    pka_precision=1.0
)
SMI = str(dimorphite_dl.protonate(SMI)[0])

mol = Chem.MolFromSmiles(SMI)
sdm = pretreat.StandardizeMol()
mol = sdm.disconnect_metals(mol)

m = Chem.MolFromSmiles(SMI)
im=Draw.MolToImage(m)

logd = scopy.ScoDruglikeness.molproperty.CalculateLogD(mol)
mr = scopy.ScoDruglikeness.molproperty.CalculateMolMR(mol)
print("logD: " + str(logd))
print("CrippenMR: " + str(mr))

tcl1 = ( ( logd - 2.009365 ) / 1.878026 ) * 1.1626256
tcl2 = ( ( mr - 90.80861 ) / 35.30108 ) * 1.9764294

tcl3 = 1 / ( 1 + ( 2.718281828459045 ** ( -1 * ( 0.7999132 + tcl1 + tcl2 ) ) ) )

print("TC/L probability: " + str(tcl3))

st.image(im)

