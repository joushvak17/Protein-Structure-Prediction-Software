from Bio import PDB


def extract_experimental(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    experimental = structure.header.get("structure_method", "Unknown")
    res = structure.header.get("resolution", "Unknown")
    return experimental, res
