# TC/L interaction probability model
Python web app and script to compute probability of interaction of a molecule with taurocholate/lecithin at pH 6.5, based on two descriptors: logD and CrippenMR. The model is inspired by Mol. Pharmaceutics 2022, 19, 2868−2876 (https://doi.org/10.1021/acs.molpharmaceut.2c00227) but was built with descriptors from rdkit/scopy instead of MOE/PaDEL, using logD for pH 7.4 instead of 7.0. For the same traning/validation set similar statistics are retrieved (balanced accuracy: 0.86/0.83, AUC: 0.93/0.93).

Statistics regarding the similarity of the compound's properties (logD/CrippenMR) and chemical structure (SDC applicability domain metric, J. Chem. Inf. Model. 2019, 59, 1, 181–189) in relation to the training set are listed to further assess predictability.

## Web application

Streamlit web app ```tclint_streamlit.py``` is accessible at:

https://tclint.streamlit.app/

## Local command-line script

To use the command-line script ```tclint.py``` locally, run the following code in the downloaded folder:
```
python tclint.py "SMILESCODE" 
```
Instead of a single SMILES code, a "*.dat" file can be given as input, with one SMILES code per line (without parentheses). In this case, the script outputs the properties for each molecule into "results.dat":
```
python tclint.py listofsmiles.dat 
```

Packages needed for running:
```
rdkit 2022.09.1
scopy 1.2.5
dimorphite-dl 1.3.2
numpy 1.23.4
```

## Conda program installation from source

The model can be installed as program in a new conda environment with all the necessary packages from source in the following way:

1. *go to your local anaconda/envs/ directory*
2. ```git clone https://github.com/juppifluppi/tclint.git```
3. ```cd tclint```
4. ```conda env create -f env.yml```
5. ```conda activate tclint```
6. ```pip install -e .```

You can then use the alias tclint:
```
tclint "SMILESCODE" 
tclint listofsmiles.dat
```

## Databases of original models for TC/L and lipid interactions

Alternatively, ```db_enamine.html``` and ```db_zinc.html``` contain already calculated TC/L and lipid interaction probabilities for drug molecules with the original models 2 and 3 from the publication, using MOE/PaDEL descriptors (not manually inspected for correct protonation!).

enamine_db.html: Enamine FDA approved drug collection (07.08.2022)

zinc_db.html: ZINC DrugBank (20.02.18)

## Please cite

[Schlauersbach, Jonas, et al. "Predicting Bile and Lipid Interaction for Drug Substances." Molecular pharmaceutics 19.8 (2022): 2868-2876.](https://doi.org/10.1021/acs.molpharmaceut.2c00227)
