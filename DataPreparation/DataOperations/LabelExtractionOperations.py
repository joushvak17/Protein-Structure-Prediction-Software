# Import the needed libraries
from Bio import PDB


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
        
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        model = structure[0]
        dssp = PDB.DSSP(model, pdb_file, dssp="/usr/local/bin/mkdssp")
        
        # Extract the secondary structure
        secondary_structure = [dssp[key][2] for key in dssp.keys()]
        
        # Extract solvent accessibility
        solvent_accessibility = [dssp[key][3] for key in dssp.keys()]

        # Extract disorder regions (assuming disorder is marked as 'X' in DSSP)
        disorder_regions = [1 if dssp[key][2] == 'X' else 0 for key in dssp.keys()]
                    
        return secondary_structure, solvent_accessibility, disorder_regions
    except ValueError:
        return None, None, None