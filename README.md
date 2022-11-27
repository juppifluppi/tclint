# TC/L interaction probability model
Python web app and script to compute probability of interaction of molecule with taurocholate/lecithin at pH 6.5, based on two descriptors: logD and CrippenMR. The model is inspired by Mol. Pharmaceutics 2022, 19, 2868−2876 (https://doi.org/10.1021/acs.molpharmaceut.2c00227) but was built with descriptors from rdkit/scopy instead of MOE/PaDEL, using logD for pH 7.4 instead of 7.0. For the same traning/validation set similar statistics are retrieved (balanced accuracy: 0.86/0.83, AUC: 0.93/0.93).

Statistics regarding the similarity of the compound's properties (logD/CrippenMR) and chemical structure (SDC applicability domain metric, J. Chem. Inf. Model. 2019, 59, 1, 181–189) in relation to the training set are listed to further assess predictability.

## Web application

Streamlit web app (tcli_streamlit.py) is accessible at:

https://tcli.streamlit.app/

## Local command-line script

To use the command-line script (tcli.py) locally, run the following code in the downloaded folder:
```
python tcli.py "SMILESCODE" 
```
Conda packages installed for running:
```
rdkit 2022.09.1
scopy 1.2.5
dimorphite-dl 1.3.2
numpy 1.23.4
pandas 1.5.1
```

## Original model databases for TC/L and lipid interactions

Alternatively, enamine_db.html and zinc_db.html contain already calculated TC/L and lipid interaction probabilities for drug molecules with the original models 2 and 3 from the publication, using MOE/PaDEL descriptors (not manually inspected for correct protonation!).
