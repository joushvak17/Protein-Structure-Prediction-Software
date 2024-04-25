from Bio import PDB


def extract_experimental(pdb_file):
    r_value = 50
    r_free = 50
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    experimental = structure.header.get("structure_method", "Unknown")
    res = structure.header.get("resolution", "Unknown")

    with open(pdb_file, "r") as f:
        for line in f:
            try:
                if "REMARK   3   R VALUE" in line:
                    r_value = float(line.split()[-1])
                elif "REMARK   3   FREE R VALUE" in line:
                    r_free = float(line.split()[-1])
            except ValueError:
                pass

    return experimental, res, r_value, r_free
