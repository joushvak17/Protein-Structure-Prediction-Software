from Bio.PDB import PDBParser
from Bio.PDB.DSSP import make_dssp_dict
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
    disorder_prediction = []
    
    for pdb_file in pdb_files:
        structure = PDBParser(QUIET=True).get_structure(pdb_file, f"DataPreparation/PDBData/{pdb_file}.pdb")
        model = structure[0]
        dssp_dict, _ = make_dssp_dict(model)
    
        disorder_prediction.append([dssp_dict[key][2] for key in dssp_dict])