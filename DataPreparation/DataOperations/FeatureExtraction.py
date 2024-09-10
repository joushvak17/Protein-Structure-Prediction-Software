# Import the needed libraries
import math

from Bio import SeqIO, AlignIO, Phylo
from Bio.Align import AlignInfo
from DataOperations.FeatureExtractionOperations import *

# NOTE: The extractions for the unaligned sequences are done
def extract_unaligned(path):
    """Function to extract the features from the unaligned sequences

    Args:
        path (str): The path to the unaligned sequences

    Returns:
        dict: The dictionary containing the extracted features
    """
    # Define the unaligned dataframe that will have all the calculated feature values
    unaligned_data = {
        "ID": [],
        "Unaligned Sequence": [],
        'A': [], 'R': [], 'N': [], 'D': [],
        'C': [], 'E': [], 'Q': [], 'G': [],
        'H': [], 'I': [], 'L': [], 'K': [],
        'M': [], 'F': [], 'P': [], 'S': [],
        'T': [], 'W': [], 'Y': [], 'V': [],
        "Hydrophobicity (Kyte-Doolittle Scale)": [],
        "Net Charge at pH 7.0 (Neutral)": [],
        "Isoelectric Point": [],
        "Molecular Weight": [],
        "Sequence Length": [],
        "Dipeptide Composition": [],
        "Helix Propensity": [],
        "Sheet Propensity": [],
        "Coil Propensity": [],
        "Disorder Prediction": []
    } 

    for seq_record in SeqIO.parse(path, "fasta"):
        # ID and Unaligned Sequence
        unaligned_data["ID"].append(seq_record.id)
        unaligned_data["Unaligned Sequence"].append(str(seq_record.seq))
        
        # Amino Acid Frequency
        aa_composition = calculate_a_acid_composition(str(seq_record.seq))
        for amino_acid, percent in aa_composition.items():
            unaligned_data[amino_acid].append(percent)
        
        # Hyrdophobicity
        hydrophobicity_values = calculate_hydrophobicity(str(seq_record.seq))
        unaligned_data["Hydrophobicity (Kyte-Doolittle Scale)"].append(hydrophobicity_values)
        
        # Net Charge at pH 7.0
        charge_7 = calculate_polarity(str(seq_record.seq), 7.0)
        unaligned_data["Net Charge at pH 7.0 (Neutral)"].append(charge_7)
        
        # Isolectric Point
        isolectric_values = calculate_pI(str(seq_record.seq))
        unaligned_data["Isoelectric Point"].append(isolectric_values)
        
        # Molecular Weight
        mw_values = calculate_mw(str(seq_record.seq))
        unaligned_data["Molecular Weight"].append(mw_values)
        
        # Sequence Length
        sequence_length = len(seq_record.seq)
        unaligned_data["Sequence Length"].append(sequence_length)
        
        # Dipeptide Composition
        dipeptide_composition = calculate_dipeptide_composition(str(seq_record.seq))
        unaligned_data["Dipeptide Composition"].append(dipeptide_composition)
        
        # Secondary Structure Propensity
        helix_propensity, sheet_propensity, coil_propensity = calculate_secondary_structure_propensity(str(seq_record.seq))
        unaligned_data["Helix Propensity"].append(helix_propensity)
        unaligned_data["Sheet Propensity"].append(sheet_propensity)
        unaligned_data["Coil Propensity"].append(coil_propensity)
        
        # Disorder Prediction
        disorder_prediction = predict_disorder(seq_record.seq, 0.5)
        unaligned_data["Disorder Prediction"].append(disorder_prediction)   
        
    return unaligned_data

# TODO: Figure out this entire function
def extract_aligned(path):
    """Function to extract the features from the aligned sequences

    Args:
        path (str): The path to the aligned sequences

    Returns:
        dict: The dictionary containing the extracted features
    """
    # Read the alignment sequences
    alignment = AlignIO.read(path, "fasta")

    # FIXME: Several of the features, such as information_content, will be deprecated
    # Calculate Conservation Score
    summary_info = AlignInfo.SummaryInfo(alignment)
    positional_conservation_scores = []
    e_freq_table = {char: 0.05 for char in "ACDEFGHIKLMNPQRSTVWY"}

    for i in range(alignment.get_alignment_length()):
        score = summary_info.information_content(start=i, end=i+1, 
                                                 e_freq_table=e_freq_table, 
                                                 chars_to_ignore=["-"])
        positional_conservation_scores.append(score)
        
    # Get the consensus sequence
    consensus = summary_info.dumb_consensus(threshold=0.7, ambiguous='X')

    # Calculate the Position Specific Score Matrix
    pssm = summary_info.pos_specific_score_matrix(chars_to_ignore=["-"])
    pssm_scores = []
    
    for score in pssm:
        pssm_scores.append(score)
    
    # Initialize variables to store gap statistics
    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)
    gap_count_per_position = [0] * alignment_length

    # Count the number of gaps at each position
    for seq_record in alignment:
        for i, residue in enumerate(str(seq_record.seq)):
            if residue == "-":
                gap_count_per_position[i] += 1
                
    # Calculate the percentage of gaps at each position
    perc_gap_per_position = [count / num_sequences * 100 for count in gap_count_per_position]

    # Calculate average gap length
    all_gaps = []
    for seq_record in alignment:
        sequence = str(seq_record.seq)
        gaps = [gap for gap in sequence.split('-') if gap]
        gaps_length = [len(gap) for gap in gaps]
        all_gaps.extend(gaps_length)
        
    # Calculate the positional entropy
    entropy_list = []
    
    for i in range(alignment_length):
        freq_dict = {}
        for seq_record in alignment:
            residue = str(seq_record.seq)[i]
            if residue in freq_dict:
                freq_dict[residue] += 1
            else:
                freq_dict[residue] = 1
                
        frequencies = [count / num_sequences for count in freq_dict.values()]
        
        entropy = -sum(f * math.log2(f) for f in frequencies if f > 0)
        entropy_list.append(entropy)
    
    # Read the phylogenetic tree
    tree_file = "DataPreparation/FASTAData/tree.newick"
    tree = Phylo.read(tree_file, "newick")
    
    # Create a dictionary to store the phylogenetic weighting
    phylo_weights = {}
    
    # Calculate the total branch length for each terminal node
    for terminal_node in tree.get_terminals():
        node_path = tree.get_path(terminal_node)
        phylo_weights[terminal_node.name] = sum(clade.branch_length for clade in node_path)
    
    # Normalize the weights
    total_weight = sum(phylo_weights.values())
    for key in phylo_weights:
        phylo_weights[key] /= total_weight

    """
    TODO: 
    Figure out how many data points are needed for each feature. Some of the features are 
    calculated per position, while others are calculated per sequence. Per sequence features
    will be about 1818, while per position features will be about 2558. Can aggregate the
    per position features to be the same length as the per sequence features.    
    """
    # Define the aligned dataframe that will have the calculated feature values
    aligned_data = {
        "ID": [], # 1818
        "Aligned Sequence": [], # 1818
        "Conservation Scores": positional_conservation_scores, # 2558
        # TODO: Figure out what the output for the consensus sequence is
        "Consensus Sequence": consensus,
        "PSSM Scores": pssm_scores, # 2558
        "Percentage of Gaps Per Position": perc_gap_per_position, # 2558
        "Positional Entropy": entropy_list, # 2558
        "Phylogenetic Weighting": phylo_weights, # 1818
        "Sequence Length": [], # 1818
        "Gap Count": [], # 1818
        "Percentage Gaps": [] # 1818
    }

    for seq_record in alignment:
        aligned_data["ID"].append(seq_record.id)
        aligned_data["Aligned Sequence"].append(str(seq_record.seq))
        
        sequence = str(seq_record.seq)
        len_sequence = len(sequence)
        gap_count = sequence.count('-')
        perc_gaps = (gap_count / len_sequence) * 100

        aligned_data["Sequence Length"].append(len_sequence)
        aligned_data["Gap Count"].append(gap_count)
        aligned_data["Percentage Gaps"].append(perc_gaps)
        
    return aligned_data