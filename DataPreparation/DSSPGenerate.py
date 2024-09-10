# Import the needed libraries
import os
import subprocess


def generate_dssp(pdb_files, pdb_data) -> None:
    """Function to generate the DSSP files

    Args:
        pdb_files (list): List of PDB files
        pdb_data (str): The path to the PDB data
    """
    # Create the DSSP directory if it does not exist
    if not os.path.exists("DataPreparation/DSSPData"):
        os.mkdir("DataPreparation/DSSPData")
    
    # Generate the DSSP files, looping through the PDB files
    for pdb_file in pdb_files:
        # Define the DSSP file and PDB file
        dssp_file = f"DataPreparation/DSSPData/{pdb_file}.dssp"
        pdb_file = f"{pdb_data}/{pdb_file}.pdb"
        
        # Generate the DSSP file
        cmd = [
            "mkdssp",
            pdb_file,
            dssp_file
        ]
        
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}, {e.output}")
        
        # Check which DSSP file was not generated
        if not os.path.exists(dssp_file):
            print(f"Failed to generate DSSP file for {pdb_file}")
            
def main() -> None:
    # Process the PDB files
    pdb_data = "DataPreparation/PDBData"
    pdb_files_with_extension = os.listdir(pdb_data)
    pdb_files = [file[:-4] for file in pdb_files_with_extension if file.endswith(".pdb")]
    
    # Generate the DSSP files
    generate_dssp(pdb_files, pdb_data)
  
if __name__ == "__main__":
    main()