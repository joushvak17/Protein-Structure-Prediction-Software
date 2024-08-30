# Import the needed libraries
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from DataOperations.LabelExtractionOperations import extract_data


def extract_labels(path):
    # Define the label dataframe that will be the prediction outputs
    # TODO: Check to see if experimental method is needed
    # - The ID and R Free value have been removed
    label_data = {"Experimental": [], 
                  "Resolution": [],
                  "R Value": []}

    with (ThreadPoolExecutor() as executor):
        for experimental, res, r_value in executor.map(extract_data, SeqIO.parse(path, "fasta")):
            label_data["Experimental"].append(experimental)
            label_data["Resolution"].append(res)
            label_data["R Value"].append(r_value)
