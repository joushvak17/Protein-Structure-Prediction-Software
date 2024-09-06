# Import the needed libraries
import requests
import os
import urllib3
import subprocess
import logging

from retrying import retry
from tqdm import tqdm
from Bio.PDB import PDBParser, is_aa, Polypeptide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

@retry(stop_max_attempt_number=5, wait_fixed=2000)
def download_pdb(pdb_id, pdb_data):
    """Function that will download the PDB files

    Args:
        pdb_id (str): The PDB ID
        pdb_data (str): The path to save the PDB files
    """
    try:
        # Define the URL to download the PDB file
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        
        # Define the path to save the PDB file
        path = pdb_data + f"/{pdb_id}.pdb"
        
        # Download the PDB file
        response = requests.get(url)
        if response.status_code == 200:
            with open(path, "wb") as f:
                f.write(response.content)
    except urllib3.exceptions.HTTPError as e:
        logging.error(f"Failed to download PDB file for {pdb_id}: {e}")
        raise

def get_pdb_id(pdb_csv_file):
    """Function to get the PDB files from the PDB IDs

    Args:
        pdb_csv_file (str): File path to the CSV file containing the PDB IDs
    """
    logging.debug("Getting PDB files from the PDB IDs")
    
    # Get the PDB IDs from the CSV file
    with open(pdb_csv_file, "r") as f:
        line = f.readline()
        pdb_ids = line.split(",")

    logging.debug(f"The total length of the IDs in the CSV file is: {len(pdb_ids)}")

    # Define a directory to save the PDB files
    pdb_data = "DataPreparation/PDBData"

    # Create the directory to save the PDB files
    if not os.path.exists(pdb_data):
        os.makedirs(pdb_data)
    
    # Download all files from PDB using the IDs that were extracted
    with tqdm(total=len(pdb_ids), desc="Downloading PDB files") as pbar:
        for pdb_id in pdb_ids:
            download_pdb(pdb_id, pdb_data)
            pbar.update(1)
    
    logging.debug("The total length of the IDs that were able to be downloaded is: %s", len(os.listdir(pdb_data)))

def preprocess_sequence(pdb_files, pdb_data):
    """Function that will preprocess the sequences

    Args:
        pdb_files (list): List of PDB files
        pdb_data (str): The path to the PDB files

    Returns:
        list: List of sequence records
    """
    logging.debug("Preprocessing the sequences")
    
    sequence_records = []
    sequences_seen = set()
    
    # Count the number of missing atoms and non-standard residues 
    non_standard_residue_count = 0
    missing_atom_count = 0
    
    # Loop through all the PDB files
    for pdb_file in pdb_files:
        
        # Parse the PDB files 
        structure = PDBParser(QUIET=True).get_structure(pdb_file, f"{pdb_data}/{pdb_file}.pdb")

        for model in structure:
            for chain in model:
                sequence = ""
                for residue in chain:

                    # Check if residue is not a standard amino acid
                    if not is_aa(residue):
                        non_standard_residue_count += 1
                        continue

                    # Check if residue is missing atoms  
                    if residue.is_disordered():
                        missing_atom_count += 1
                        continue

                    # Convert to one-letter code
                    # NOTE: There is no 'three_to_one' method in the Polypeptide module
                    # so we use three_to_index and index_to_one instead
                    # "AttributeError: module 'Bio.PDB.Polypeptide' has no attribute 'three_to_one'. 
                    # Did you mean: 'three_to_index'?"
                    try:
                        sequence += Polypeptide.index_to_one(Polypeptide.three_to_index(residue.get_resname()))
                    except KeyError:
                        continue
                        
                # Add correctly formed sequences, check if sequence is not empty and not seen before
                if sequence and sequence not in sequences_seen:
                    sequence_records.append(SeqRecord(Seq(sequence), id=f"{pdb_file}_{chain.id}",
                                                      description=f"Source File: {pdb_file}, Chain: {chain.id}"))
                    sequences_seen.add(sequence)

    logging.debug((f"Found {non_standard_residue_count} non-standard residues"))
    logging.debug((f"Found {missing_atom_count} missing atoms"))

    return sequence_records

def main():
    logging.debug("Starting main function")
    
    # Define the path to the csv file containing the PDB IDs
    pdb_csv_file = "DataPreparation/PDBDataID.csv"
    
    # Get the PDB files from the PDB IDs
    get_pdb_id(pdb_csv_file)
    
    # Preprocess the sequences
    pdb_data = "DataPreparation/PDBData"
    pdb_files_with_extension = os.listdir(pdb_data)
    pdb_files = [file[:-4] for file in pdb_files_with_extension if file.endswith(".pdb")] 
    sequences = preprocess_sequence(pdb_files, pdb_data)
    
    # Define a path to save the FASTA file
    fasta_path = "DataPreparation/FASTAData"
    
    # Create the directory to save the FASTA file
    if not os.path.exists(fasta_path):
        os.makedirs(fasta_path)

    # Write the sequences to a FASTA file
    SeqIO.write(sequences, fasta_path + "/Sequences.fasta", "fasta")

    # Align the sequences using Clustal Omega
    in_file = fasta_path + "/Sequences.fasta"
    out_file = fasta_path + "/Aligned_Sequences.fasta"
    
    if not os.path.exists(in_file):
        logging.error(f"Error: Input file {in_file} does not exist.")
    elif os.path.getsize(in_file) == 0:
        logging.error(f"Error: Input file {in_file} is empty.")
    else:
        logging.debug(f"Aligning sequences using Clustal Omega")
        
        # Define the command to run Clustal Omega
        cmd = [
            "clustalo",
            "-i", in_file,
            "-o", out_file,
            "--auto",
            "--threads", str(os.cpu_count()),
            "--verbose"
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error: {e}, {e.output}")

if __name__ == "__main__":
    main()