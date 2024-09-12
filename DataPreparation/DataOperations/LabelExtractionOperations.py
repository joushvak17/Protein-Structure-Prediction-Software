# Import the needed libraries
import os

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def extract_data(seq_record):
    """Function to extract the labels from the PDB file

    Args:
        seq_record (str): The sequence record

    Returns:
        tuple: The secondary structure, solvent accessibility, and disorder regions
    """
    try:
        seq_id = seq_record.id
        base_pdb_id = seq_id.split("_")[0]
        pdb_file = f"DataPreparation/PDBData/{base_pdb_id}.pdb"
        dssp_file = f"DataPreparation/DSSPData/{base_pdb_id}.dssp"
        
        if os.path.exists(dssp_file): 
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(base_pdb_id, pdb_file)
            model = structure[0]
            dssp = DSSP(model, dssp_file, dssp="mkdssp", file_type="DSSP")
            
            # Extract the secondary structure
            secondary_structure = [dssp[key][2] for key in dssp.keys()]
            
            # Extract solvent accessibility
            solvent_accessibility = [dssp[key][3] for key in dssp.keys()]

            # Extract disorder regions (assuming disorder is marked as 'X' in DSSP)
            disorder_regions = [1 if dssp[key][2] == 'X' else 0 for key in dssp.keys()]
                        
            return secondary_structure, solvent_accessibility, disorder_regions
        else:
            return None, None, None
    except ValueError:
        return None, None, None