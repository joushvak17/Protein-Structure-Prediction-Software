from Bio import PDB

def extract_data(seq_record):
        try:
            r_value = None
            seq_id = seq_record.id
            base_pdb_id = seq_id.split("_")[0]
            pdb_file = f"DataPreparation/PDBData/{base_pdb_id}.pdb"
            
            try:
                parser = PDB.PDBParser(QUIET=True)
                structure = parser.get_structure("protein", pdb_file)
                experimental = structure.header.get("structure_method", "Unknown")
                res = structure.header.get("resolution", "Unknown")

            except ValueError:
                experimental, res = None, None
            
            with open(pdb_file, "r") as f:
                for line in f:
                    try:
                        if "REMARK   3   R VALUE            (WORKING SET) :" in line:
                            r_value = float(line.split()[-1])
                    except ValueError:
                        pass
                        
            return experimental, res, r_value
        except ValueError:
            return None, None, None