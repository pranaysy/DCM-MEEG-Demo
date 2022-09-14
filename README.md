# Tutorial for group Dynamic Causal Modelling (DCM) of M/EEG
Demo scripts for Dynamic Causal Modelling (DCM) of M/EEG data from the Wakeman & Henson (2015) open dataset using a hierarchical Bayesian framework called Parametric Empirical Bayes (PEB). 

### Software
- Clone this repository on your system. 
- The tutorial uses the current stable version of SPM12 with minor patches available here: https://github.com/pranaysy/spm12  Please set this as the primary SPM version on your MATLAB path for use with this tutorial.

### Dataset
The raw BIDS-formatted dataset can be obtained from [OpenNeuro](https://openneuro.org/datasets/ds000117/versions/1.0.5) using tools like [datalad](https://www.datalad.org/) or git-annex. Processing of the raw data is done according to Henson et al (2019), and a [convenience script is provided](https://github.com/pranaysy/DCM-MEEG-Demo/blob/094c28fa49dbaa21fc1b41451e74cd6326c7c30f/code/spm_master_script_data_preprocessing.m) for processing data in preparation for this tutorial. Additionally, tutorial-ready data with forward models can be downloaded from figshare.

Once processed, the data should be stored as indicated in [./data/folder_structure.txt](https://github.com/pranaysy/DCM-MEEG-Demo/blob/094c28fa49dbaa21fc1b41451e74cd6326c7c30f/data/folder_structure.txt)

### References
1. Wakeman, D.G. & Henson, R.N. (2015). A multi-subject, multi-modal human neuroimaging dataset. Sci. Data 2:150001 doi: 10.1038/sdata.2015.1
2. Henson RN, Abdulrahman H, Flandin G and Litvak V (2019) Multimodal Integration of M/EEG and f/MRI Data in SPM12. Front. Neurosci. 13:300. doi: 10.3389/fnins.2019.00300
