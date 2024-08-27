# Import the needed libraries
from collections import Counter
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

# TODO: Figure out how this function works
def calculate_dipeptide_composition(sequence):
    dipeptides = [sequence[i:i+2] for i in range(len(sequence) - 1)]
    dipeptide_counts = Counter(dipeptides)
    total_dipeptides = sum(dipeptide_counts.values())
    dipeptide_freq = {dipeptide: count / total_dipeptides for dipeptide, count in dipeptide_counts.items()}
    return dipeptide_freq

def predict_disorder(sequence, threshold):
    disorder_scores = meta.predict_disorder(sequence)
    average_disorder_score = sum(disorder_scores) / len(disorder_scores)
    is_disordered = 1 if average_disorder_score > threshold else 0
    return is_disordered