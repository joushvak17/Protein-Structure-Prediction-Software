from Bio import PDB


def extract_experimental(pdb_file):
    parser = PDB.PDBParser()
    try:
        structure = parser.get_structure("protein", pdb_file)
        experimental = structure.header.get("structure_method", "Unknown")
    except Exception as e:
        print(f"Error parsing {pdb_file}: {str(e)}")
        experimental = "Error"
    return experimental
