from Bio.PDB import PDBParser, DSSP
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def calculate_a_acid_composition(sequence):
    protein_analysis = ProteinAnalysis(sequence)
    aac_dict = protein_analysis.get_amino_acids_percent()
    return aac_dict


def calculate_hydrophobicity(sequence):
    hydrophobicity_values = ProteinAnalysis(sequence).gravy()
    return hydrophobicity_values


def calculate_polarity(sequence, ph):
    protein_analysis = ProteinAnalysis(sequence)
    net_charge = protein_analysis.charge_at_pH(ph)
    return net_charge


def calculate_mw(sequence):
    protein_analysis = ProteinAnalysis(sequence)
    molecular_weight = protein_analysis.molecular_weight()
    return molecular_weight


def calculate_pI(sequence):
    protein_analysis = ProteinAnalysis(sequence)
    pI = protein_analysis.isoelectric_point()
    return pI

def disorder_prediction(pdb_files):
    disorder_predictions = []
    
    for pdb_file in pdb_files:
        # Parse PDB file and create DSSP object
        structure = PDBParser(QUIET=True).get_structure(pdb_file, f"DataPreparation/PDBData/{pdb_file}.pdb")
        model = structure[0]
        dssp = DSSP(model, f"DataPreparation/PDBData/{pdb_file}.pdb")
        
        # Interpret DSSP output
        for key in dssp:
            ss, acc = dssp[key][2], dssp[key][3]
            # Here we assume coil ('C') structures with high accessibility might be disordered
            if ss == 'C' and acc > 50:  # Threshold for accessibility
                disorder_predictions.append(1)
            else:
                disorder_predictions.append(0)