# Import the needed libraries
import pandas as pd
import pickle
import os
import subprocess

from sklearn.preprocessing import LabelEncoder, StandardScaler
from DataOperations.FeatureExtraction import *
from DataOperations.LabelExtraction import *


# Define the function to save the state
def save_state(state, filename):
    with open(filename, "wb") as f:
        pickle.dump(state, f)

# Define the function to load the state
def load_state(filename):
    with open(filename, "rb") as f:
        return pickle.load(f)

def main():
    if os.path.exists(TREE_FILE) and os.path.exists(STATE_FILE):
        print("The tree file and state file exist.")
        
        # Load the saved state
        state = load_state(STATE_FILE)
        unaligned_data = state["unaligned_data"]
    else:
        # Define the path for the aligned sequences
        aligned_path = "DataPreparation/FASTAData/Aligned_Sequences.fasta"
        
        # Define the command to run Clustal Omega
        cmd = ["clustalo", "-i", aligned_path, "--guidetree-out", TREE_FILE]
            
        try:
            # Run the command
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}, {e.output}")
        
        # Define the path for the unaligned sequences
        unaligned_path = "DataPreparation/FASTAData/Sequences.fasta"
    
        # Extract the unaligned data
        unaligned_data = extract_unaligned(unaligned_path)
        
        # Save the state
        state = {"unaligned_data": unaligned_data}
        save_state(state, STATE_FILE)
        print("State saved. Please run the script again.")
        return
    
    # Define the path for the unaligned and aligned sequences
    unaligned_path = "DataPreparation/FASTAData/Sequences.fasta"
    aligned_path = "DataPreparation/FASTAData/Aligned_Sequences.fasta"
        
    aligned_data = extract_aligned(aligned_path)
    label_data = extract_labels(unaligned_path)

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
    
    # Save the final state
    state= {"aligned_data": aligned_data, 
            "label_data": label_data,
            "df": df,
            "le": le,
            "scaler": scaler}

if __name__ == "__main__":
    # Define the tree file
    TREE_FILE = "DataPreparation/FASTAData/tree.newick"

    # Define the state file
    STATE_FILE = "DataPreparation/state.pkl"
    
    main()