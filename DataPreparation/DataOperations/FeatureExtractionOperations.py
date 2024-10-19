# Import the needed libraries
# TODO: Import the metapredict library when it is available
# import metapredict as meta

from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def calculate_a_acid_composition(sequence):
    """Function to calculate the amino acid composition of a sequence

    Args:
        sequence (str): The sequence to calculate the amino acid composition

    Returns:
        dict: The dictionary containing the amino acid composition
    """
    protein_analysis = ProteinAnalysis(sequence)
    aac_dict = protein_analysis.get_amino_acids_percent()
    return aac_dict

def calculate_hydrophobicity(sequence):
    """Function to calculate the hydrophobicity of a sequence

    Args:
        sequence (str): The sequence to calculate the hydrophobicity

    Returns:
        float: The hydrophobicity value
    """
    hydrophobicity_values = ProteinAnalysis(sequence).gravy()
    return hydrophobicity_values

def calculate_polarity(sequence, ph):
    """Function to calculate the polarity of a sequence at a given pH

    Args:
        sequence (str): The sequence to calculate the polarity
        ph (float): The pH value to calculate the polarity, usually 7.0

    Returns:
        float: The polarity value
    """
    protein_analysis = ProteinAnalysis(sequence)
    net_charge = protein_analysis.charge_at_pH(ph)
    return net_charge

def calculate_mw(sequence):
    """Function to calculate the molecular weight of a sequence

    Args:
        sequence (str): The sequence to calculate the molecular weight

    Returns:
        float: The molecular weight value
    """
    protein_analysis = ProteinAnalysis(sequence)
    molecular_weight = protein_analysis.molecular_weight()
    return molecular_weight

def calculate_pI(sequence):
    """Function to calculate the isoelectric point of a sequence

    Args:
        sequence (str): The sequence to calculate the isoelectric point

    Returns:
        float: The isoelectric point value
    """
    protein_analysis = ProteinAnalysis(sequence)
    pI = protein_analysis.isoelectric_point()
    return pI

# TODO: Figure out how this function works
def calculate_dipeptide_composition(sequence):
    """Function to calculate the dipeptide composition of a sequence

    Args:
        sequence (str): The sequence to calculate the dipeptide composition

    Returns:
        dict: The dictionary containing the dipeptide composition
    """
    dipeptides = [sequence[i:i+2] for i in range(len(sequence) - 1)]
    dipeptide_counts = Counter(dipeptides)
    total_dipeptides = sum(dipeptide_counts.values())
    dipeptide_freq = {dipeptide: count / total_dipeptides for dipeptide, count in dipeptide_counts.items()}
    return dipeptide_freq

# FIXME: The values in the propensities dictionary should be replaced with actual values
def calculate_secondary_structure_propensity(sequence):
    """Function to calculate the secondary structure propensity of a sequence

    Args:
        sequence (str): The sequence to calculate the secondary structure propensity

    Returns:
        tuple: The tuple containing the helix propensity, sheet propensity, and coil propensity
    """
    # Example propensities (values are illustrative)
    propensities = {
        'H': {'A': 1.45, 'R': 0.97, 'N': 0.67, 'D': 1.01, 'C': 0.77, 'E': 1.53, 'Q': 1.23, 'G': 0.57, 'H': 1.00, 'I': 1.08, 'L': 1.34, 'K': 1.07, 'M': 1.20, 'F': 1.12, 'P': 0.59, 'S': 0.79, 'T': 0.82, 'W': 1.14, 'Y': 0.61, 'V': 1.06},
        'E': {'A': 0.97, 'R': 0.93, 'N': 0.65, 'D': 0.80, 'C': 1.30, 'E': 0.26, 'Q': 0.90, 'G': 0.81, 'H': 0.87, 'I': 1.60, 'L': 1.22, 'K': 0.74, 'M': 1.67, 'F': 1.28, 'P': 0.62, 'S': 0.72, 'T': 1.20, 'W': 1.19, 'Y': 1.29, 'V': 1.65},
        'C': {'A': 0.77, 'R': 0.90, 'N': 1.56, 'D': 0.72, 'C': 1.19, 'E': 0.47, 'Q': 0.70, 'G': 1.64, 'H': 0.87, 'I': 0.47, 'L': 0.53, 'K': 0.74, 'M': 0.60, 'F': 0.59, 'P': 1.52, 'S': 1.56, 'T': 1.20, 'W': 0.54, 'Y': 1.29, 'V': 0.79}
    }
    helix_propensity = sum(propensities['H'].get(aa, 0) for aa in sequence) / len(sequence)
    sheet_propensity = sum(propensities['E'].get(aa, 0) for aa in sequence) / len(sequence)
    coil_propensity = sum(propensities['C'].get(aa, 0) for aa in sequence) / len(sequence)
    return helix_propensity, sheet_propensity, coil_propensity

# # TODO: Figure out if this can be replaced with another tool
# def predict_disorder(sequence, threshold):
#     """Function to predict if a sequence is disordered

#     Args:
#         sequence (str): The sequence to predict if it is disordered
#         threshold (float): The threshold to determine if a sequence is disordered, usually 0.5

#     Returns:
#         Literal[0, 1]: The value to determine if a sequence is disordered, 1 if disordered, 0 if not
#     """
#     disorder_scores = meta.predict_disorder(sequence)
#     average_disorder_score = sum(disorder_scores) / len(disorder_scores)
#     is_disordered = 1 if average_disorder_score > threshold else 0
#     return is_disordered