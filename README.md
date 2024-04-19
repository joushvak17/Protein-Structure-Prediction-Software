# Protein-Structure-Prediction-Software

Protein structure prediction software that allows users to visualize predicted protein structures through the use of BioPython, OpenMM, and PyTorch.

## File Outline

Data Preparation:

- DataCollection.ipynb: File that performs the data collection of PDB files from the PDBDataID.csv file
- DataPreprocessing.ipynb: File that performs the data preprocessing of PDB files, which includes retrieval, filtering/cleaning, and alignment
- FeatureExtraction.ipynb: File that performs the feature extractions for the aligned and unaligned sequences
- FeatureExtractionOperations.py: File that has the defined functions for the feature extraction
- PDBData Folder: Contains the files that has the PDB data that this project will perform training/analysis on
- PDBDataID.csv: Contains the ID's of the PDB data that this project will perform training/analysis on
- aligned_sequences.fasta: FASTA file that has the sequences in an aligned format
- sequences.fasta: FASTA file that has the sequences in an unaligned format
