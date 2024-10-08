# Import the needed libraries
import pandas as pd
import pickle
import os
import subprocess
import logging

from DataOperations.FeatureExtraction import *
from DataOperations.LabelExtraction import *

# Configure logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def save_state(state, filename):
    """Function to save the state

    Args:
        state (dict): The state to save
        filename (str): The filename to save the state
    """
    with open(filename, "wb") as f:
        pickle.dump(state, f)

def load_state(filename):
    """Function to load the state

    Args:
        filename (str): The filename to load the state

    Returns:
        state: The loaded state
    """
    with open(filename, "rb") as f:
        return pickle.load(f)

def main():
    # Define the tree file
    TREE_FILE = "DataPreparation/FASTAData/tree.newick"

    # Define the state file
    STATE_FILE = "DataPreparation/state.pkl"
    
    if os.path.exists(TREE_FILE) and os.path.exists(STATE_FILE):
        logging.debug("The tree file and state file exist.")
        
        # Load the saved state
        state = load_state(STATE_FILE)
        unaligned_data = state["unaligned_data"]
        aligned_data = state["aligned_data"]
        label_data = state["label_data"]
    else:
        logging.debug("Tree file or state file does not exist. Running FastTree.")
        
        # Define the path for the unaligned and aligned sequences
        unaligned_path = "DataPreparation/FASTAData/Sequences.fasta"
        aligned_path = "DataPreparation/FASTAData/Aligned_Sequences.fasta"
        
        # Define the command to run FastTree
        cmd = [
            "FastTree", 
            aligned_path, 
        ]
            
        try:
            # Run the command
            with open(TREE_FILE, "w") as f:
                subprocess.run(cmd, stdout=f, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error: {e}, {e.output}")
            return
        
        logging.debug("FastTree completed. Will now extract the sequences.")
    
        # Extract the data
        unaligned_data = extract_unaligned(unaligned_path)
        aligned_data = extract_aligned(aligned_path)
        label_data = extract_labels(unaligned_path)
        
        # Save the state
        state = {
            "unaligned_data": unaligned_data,
            "aligned_data": aligned_data,
            "label_data": label_data
        }
        save_state(state, STATE_FILE)
        logging.debug("State saved. Please run the script again.")
        return
    
    # Create a save folder
    if not os.path.exists("DataPreparation/CSVData"):
        os.makedirs("DataPreparation/CSVData")
    
    # Convert the dictionaries to dataframes and save them to csv files
    unaligned_df = pd.DataFrame(unaligned_data)
    unaligned_df.to_csv("DataPreparation/CSVData/UnalignedData.csv", index=False)
    
    aligned_df = pd.DataFrame(aligned_data)
    aligned_df.to_csv("DataPreparation/CSVData/AlignedData.csv", index=False)
    
    label_df = pd.DataFrame(label_data)
    label_df.to_csv("DataPreparation/CSVData/LabelData.csv", index=False)

if __name__ == "__main__":
    main()