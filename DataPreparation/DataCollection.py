# Import the needed libraries
import requests
import os
import urllib3
import subprocess

from retrying import retry
from Bio.PDB import PDBParser, is_aa, Polypeptide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# Define a function that will download a PDB file given a PDB ID
@retry(stop_max_attempt_number=5, wait_fixed=2000)
def download_pdb(pdb_id):
    try:
        # Define the URL to download the PDB file
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        
        # Create the directory to save the PDB files
        if not os.path.exists("DataPreparation/PDBData"):
            os.makedirs("DataPreparation/PDBData")
        
        # Define the path to save the PDB file
        path = f"DataPreparation/PDBData/{pdb_id}.pdb"
        
        # Download the PDB file
        response = requests.get(url)
        if response.status_code == 200:
            with open(path, "wb") as f:
                f.write(response.content)
    except urllib3.exceptions.HTTPError as e:
        print(f"Failed to download PDB file for {pdb_id}: {e}")
        raise

def get_pdb_id():
    # Get the PDB IDs from the CSV file
    with open("DataPreparation/PDBDataID.csv", "r") as f:
        line = f.readline()
        pdb_ids = line.split(",")

    # Print the total number of PDB IDs
    print("The total length of the IDs in the PDBDataID.csv file is: ", len(pdb_ids))

    # Download all files from PDB using the IDs that were extracted
    for pdb_id in pdb_ids:
        download_pdb(pdb_id)
    
    print("The total length of the IDs that were able to be downloaded is: ", len(os.listdir("DataPreparation/PDBData")))

# Define a function that will preprocess the sequences
def preprocess_sequence(pdb_files):
    sequence_records = []
    sequences_seen = set()
    
    # Count the number of missing atoms and non-standard residues 
    non_standard_residue_count = 0
    missing_atom_count = 0
    
    # Loop through all the PDB files
    for pdb_file in pdb_files:
        
        # Parse the PDB files 
        structure = PDBParser(QUIET=True).get_structure(pdb_file, f"DataPreparation/PDBData/{pdb_file}.pdb")

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
                        
                # Add correctly formed sequences
                if sequence and sequence not in sequences_seen:
                    sequence_records.append(SeqRecord(Seq(sequence), id=f"{pdb_file}_{chain.id}",
                                                      description=f"Source File: {pdb_file}, Chain: {chain.id}"))
                    sequences_seen.add(sequence)

    print(f"Found {non_standard_residue_count} non-standard residues")
    print(f"Found {missing_atom_count} missing atoms")

    return sequence_records

def main():
    get_pdb_id()
    
    # Preprocess the sequences
    pdb_files_with_extension = os.listdir("DataPreparation/PDBData")
    pdb_files = [file[:-4] for file in pdb_files_with_extension if file.endswith(".pdb")] 
    sequences = preprocess_sequence(pdb_files)
    
    # Create the directory to save the FASTA file
    if not os.path.exists("DataPreparation/FASTAData"):
        os.makedirs("DataPreparation/FASTAData")

    # Write the sequences to a FASTA file
    SeqIO.write(sequences, "DataPreparation/FASTAData/Sequences.fasta", "fasta")

    # Align the sequences using Clustal Omega
    in_file = "DataPreparation/FASTAData/Sequences.fasta"
    out_file = "DataPreparation/FASTAData/Aligned_Sequences.fasta"
    
    if not os.path.exists(in_file):
        print(f"Error: Input file {in_file} does not exist.")
    elif os.path.getsize(in_file) == 0:
        print(f"Error: Input file {in_file} is empty.")
    else:
        # Define the command to run Clustal Omega
        cmd = [
            "clustalo",
            "-i", in_file,
            "-o", out_file,
            "--auto",
            "--threads", str(os.cpu_count()),
            "--verbose",
        ]

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}, {e.output}")

if __name__ == "__main__":
    main()