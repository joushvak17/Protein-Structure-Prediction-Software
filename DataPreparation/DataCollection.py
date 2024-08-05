# Import the needed libraries
import requests
import os
import urllib3
import subprocess

from retrying import retry
from Bio.PDB import PDBParser, is_aa, Polypeptide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO

# Define a function that will download a PDB file given a PDB ID
def download_pdb(pdb_id):
    # Define the URL to download the PDB file
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    # Define the path to save the PDB file
    path = f"PDBData/{pdb_id}.pdb"
    
    # Download the PDB file
    response = requests.get(url)
    if response.status_code == 200:
        with open(path, "wb") as f:
            f.write(response.content)
    else:
        print(f"Failed to download PDB file for {pdb_id}")

# Define a function that will implement a retry strategy up to 5 times
@retry(stop_max_attempt_number=5, wait_fixed=2000)
def download_pdb_retry(pdb_id):
    try:
        download_pdb(pdb_id)
    except urllib3.exceptions.HTTPError as e:
        print(f"Failed to download PDB file for {pdb_id}: {e}")
        raise

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
        structure = PDBParser(QUIET=True).get_structure(pdb_file, f"PDBData/{pdb_file}.pdb")

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
                    try:
                        sequence += Polypeptide.three_to_one(residue.get_resname())
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
    # Get the PDB IDs from the CSV file
    with open("PDBDataID.csv", "r") as f:
        line = f.readline()
        pdb_ids = line.split(",")

    # Print the total number of PDB IDs
    print("The total length of the IDs in the PDBDataID.csv file is: ", len(pdb_ids))

    # Download all files from PDB using the IDs that were extracted
    for pdb_id in pdb_ids:
        download_pdb_retry(pdb_id)
    
    print("The total length of the IDs that were able to be downloaded is: ", len(os.listdir("PDBData")))

    # Preprocess the sequences
    pdb_files_with_extension = os.listdir("PDBData")
    pdb_files = [file[:-4] for file in pdb_files_with_extension if file.endswith(".pdb")] 
    sequences = preprocess_sequence(pdb_files)

    # Write the sequences to a FASTA file
    SeqIO.write(sequences, "FASTAData/Sequences.fasta", "fasta")

    # Align the sequences using Clustal Omega
    in_file = "FASTAData/Sequences.fasta"
    out_file = "FASTAData/Aligned_Sequences.fasta"
    clustal_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)

    try:
        subprocess.run(str(clustal_cline), check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}, {e.output}")

if __name__ == "__main__":
    main()