# Import the needed libraries
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def main() -> None:
    try:
        dssp_file = "DataPreparation/DSSPData/3G6E.dssp"
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("3G6E", "DataPreparation/PDBData/3G6E.pdb")
        model = structure[0]
        dssp = DSSP(model, dssp_file, dssp="mkdssp", file_type="DSSP")
        
        key = list(dssp.keys())[2]
        print(f"DSSP keys for 3G6E: {key}")
    except ValueError:
        print("An error occurred")
        
if __name__ == "__main__":
    main()