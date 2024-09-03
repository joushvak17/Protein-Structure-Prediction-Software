# Import the needed libraries
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from DataOperations.LabelExtractionOperations import extract_data


def extract_labels(path):
    """Function to extract the labels from the FASTA file

    Args:
        path (str): The path to the FASTA file
    """
    # Define the label dataframe that will be the prediction outputs
    label_data = {"Secondary Structure": [],
                  "Solvent Accessibility": [],
                  "Disorder Regions": []}
    
    seq_record = SeqIO.read(path, "fasta")
    
    with (ThreadPoolExecutor() as executor):
        for secondary_structure, solvent_accessibility, disorder_regions in executor.map(extract_data, seq_record):
            label_data["Secondary Structure"].append(secondary_structure)
            label_data["Solvent Accessibility"].append(solvent_accessibility)
            label_data["Disorder Regions"].append(disorder_regions)
