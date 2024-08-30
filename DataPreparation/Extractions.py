# Import the needed libraries
import pandas as pd

from sklearn.preprocessing import LabelEncoder, StandardScaler
from DataOperations.FeatureExtraction import *
from DataOperations.LabelExtraction import *


def main():
    # Define the path for the unaligned and aligned sequences
    unaligned_path = "DataPreparation/FASTAData/Sequences.fasta"
    aligned_path = "DataPreparation/FASTAData/Aligned_Sequences.fasta"
    
    # Extract the unaligned, aligned, and label data
    unaligned_data = extract_unaligned(unaligned_path)
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

if __name__ == "__main__":
    main()