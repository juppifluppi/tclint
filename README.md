# TC/L interaction probability model
Python web app and script to compute probability of interaction of molecule with taurocholate/lecithin at pH 6.5, based on two descriptors: logD and CrippenMR. The model is inspired by Mol. Pharmaceutics 2022, 19, 2868−2876 (https://doi.org/10.1021/acs.molpharmaceut.2c00227) but was built with descriptors from rdkit/scopy instead of MOE/PaDEL, using logD for pH 7.4 instead of 7.0. For the same traning/validation set similar statistics are retrieved (balanced accuracy: 0.86/0.83, AUC: 0.93/0.93).

Streamlit web app: https://tcl-interaction2.streamlit.app/

Python command-line script for download:

Lipid interaction probability model (model3.py) is currently work in progress.

Alternatively, enamine_db.html and zinc_db.html contain already calculated TC/L and lipid interaction probabilities for drug molecules with the original models 2 and 3 from the publication, using MOE/PaDEL descriptors (not manually inspected for correct protonation!).
