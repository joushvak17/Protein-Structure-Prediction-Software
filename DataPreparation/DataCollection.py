# Import the needed libraries
import requests
import os
import urllib3
from retrying import retry

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

# Get the PDB IDs from the CSV file
with open("PDBDataID.csv", "r") as f:
    line = f.readline()
    pdb_ids = line.split(",")

# Print the total number of PDB IDs
print("The total length of the IDs in the PDBDataID.csv file is: ", len(pdb_ids))

# Download all files from PDB using the IDs that were extracted
for pdb_id in pdb_ids:
    download_pdb_retry(pdb_id)
    
print("The total length of the IDs that were able to be downloaded is: ", len(os.listdir("PDB Data")))