# Import the needed libraries
import metapredict as meta

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

def predict_disorder(sequence):
    disorder_prediction = meta.predict_disorder(sequence)
    return disorder_prediction