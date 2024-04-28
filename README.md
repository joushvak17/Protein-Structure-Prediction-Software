# Protein Structure Prediction Software

Protein structure prediction software that allows users to visualize predicted protein structures through the use of 
BioPython, OpenMM, PyMol, and PyTorch.

## File Outline

Data Preparation:
- DataCollection.ipynb: File that performs the data collection with the given PDB IDs
- DataPreprocessing.ipynb: File that performs the data preprocessing and cleaning with the given PDB files
- DataSplitting.ipynb: File that performs the data splitting and does some conversion types for the columns
- FeatureExtraction.ipynb: File that performs the feature and label extractions with the given data
- LabelExtractionTest.ipynb: File that contains some tests for implementing label extraction
---
- FeatureExtractionOperations.py: File that contains the functions used to perform the feature extractions
- LabelExtractionOperations.py: File that contains the functions used to perform the label extractions
---
- Dataset.csv: File that has the extracted dataset that will be used for machine learning
- PDBDataID.csv: File that has the PDB IDs used for data collection
--- 
- sequences.fasta: File that has the unaligned sequences which are preprocessed
- aligned_sequences.fasta: File that has the aligned sequences which were done using Clustal Omega

IntegrationOpenMM:
- IntegrationOpenMM.ipynb: File that integrates OpenMM into the protein structure prediction software
---
- PDBTest.pdb: File that contains the PDB test file used for integration testing
