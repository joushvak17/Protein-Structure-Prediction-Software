# Import the needed libraries
import pandas as pd

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from concurrent.futures import ThreadPoolExecutor
from sklearn.preprocessing import LabelEncoder, StandardScaler
from DataOperations.LabelExtractionOperations import extract_experimental
from DataOperations.FeatureExtractionOperations import *

def extract_unaligned(path):
    # Define the unaligned dataframe that will have all the calculated feature values
    # - The Net Charge at pH 3.0 and 11.0 have beeen removed
    unaligned_data = {"ID": [], 
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
                    "Sequence Length": []} 

    for seq_record in SeqIO.parse(path, "fasta"):
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
        
    return unaligned_data

def extract_data(seq_record):
        try:
            r_value = None
            seq_id = seq_record.id
            base_pdb_id = seq_id.split("_")[0]
            pdb_file = f"DataPreparation/PDBData/{base_pdb_id}.pdb"
            experimental, res = extract_experimental(pdb_file)
            
            with open(pdb_file, "r") as f:
                for line in f:
                    try:
                        if "REMARK   3   R VALUE            (WORKING SET) :" in line:
                            r_value = float(line.split()[-1])
                    except ValueError:
                        pass
                        
            return experimental, res, r_value
        except ValueError:
            return None, None, None, None

def main():
    # Define the path for the unaligned sequences
    path = "DataPreparation/FASTAData/Sequences.fasta"
    
    # Extract the unaligned data
    unaligned_data = extract_unaligned(path)

    # Read the alignment sequences
    alignment = AlignIO.read("DataPreparation/FASTAData/Aligned_Sequences.fasta", "fasta")

    # Calculate Conservation Score
    start = 0
    end = alignment.get_alignment_length()
    e_freq_table = {char: 0.05 for char in "ACDEFGHIKLMNPQRSTVWY"}
    conservation_score = AlignInfo.SummaryInfo(alignment).information_content(start, end, 
                                                                              e_freq_table=e_freq_table, 
                                                                              chars_to_ignore=["-"])

    # Initialize variables to store gap statistics
    alignment_length = end
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

    # Define the aligned dataframe that will have the calculated feature values
    # TODO: Check to see if the following features are needed
    # - The Consensus Sequence, total gaps in alignment, average gap length, and mutations from consensus have been removed
    aligned_data = {"ID": [], 
                    "Aligned Sequence": [],
                    "Conservation Scores": [conservation_score] * num_sequences,
                    "Percentage of Gaps Per Position": [perc_gap_per_position] * num_sequences,
                    "Sequence Length": [],
                    "Gap Count": [],
                    "Percentage Gaps": []}

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

    # Define the label dataframe that will be the prediction outputs
    # TODO: Check to see if experimental method is needed
    # - The ID and R Free value have been removed
    label_data = {"Experimental": [], 
                  "Resolution": [],
                  "R Value": [],}

    with (ThreadPoolExecutor() as executor):
        for experimental, res, r_value in executor.map(extract_data, SeqIO.parse("DataPreparation/FASTAData/Sequences.fasta", "fasta")):
            label_data["Experimental"].append(experimental)
            label_data["Resolution"].append(res)
            label_data["R Value"].append(r_value)

    # Convert the dictionaries to dataframes
    unaligned_df = pd.DataFrame(unaligned_data)
    aligned_df = pd.DataFrame(aligned_data)
    label_df = pd.DataFrame(label_data)

    # Concatenate the dataframes
    df = pd.concat([unaligned_df, aligned_df, label_df], axis=1)
    
    # Remove the rows with missing values
    df.dropna(inplace=True)

    # Encode the experimental methods and normalize the data
    le = LabelEncoder()
    df["Experimental"] = le.fit_transform(df["Experimental"])

    scaler = StandardScaler()
    num_cols = df.select_dtypes(include=["int64", "float64"]).columns
    df[num_cols] = scaler.fit_transform(df[num_cols])

    # Print the length of the dataframe
    print(f"Length of the dataframe: {len(df)}")

    # Save the dataframe to a csv file
    df.to_csv("DataPreparation/Features.csv", index=False)

if __name__ == "__main__":
    main()