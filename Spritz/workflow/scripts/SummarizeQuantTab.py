import sys
import os
import numpy as np
import pandas as pd

def read_tpm_file(file_path):
    """Reads a TPM file and returns the Gene IDs and TPMs."""
    table = pd.read_csv(file_path, sep="\t")
    return table["Gene ID"].values, table["TPM"].values

def main():
    if len(sys.argv) < 3:
        print("Usage: SummarizeQuantTab.py <output_file> <input_file_1> <input_file_2> ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    ids_list = []
    tpms_list = []

    for file in input_files:
        currIds, currTpms = read_tpm_file(file)
        
        if not ids_list:
            ids_list = currIds
        elif not np.array_equal(currIds, ids_list):
            print(f"Error with IDs in file: {os.path.basename(file)}")
            sys.exit(1)

        tpms_list.append(currTpms)

    tpms_list.insert(0, ids_list)
    dataframe = np.row_stack(tpms_list)
    dataframe[1:, 0] = [os.path.basename(file).split(".")[0] for file in input_files]
    
    # Create a DataFrame and save to CSV
    pddf = pd.DataFrame(dataframe[1:, 1:], index=dataframe[1:, 0], columns=dataframe[0, 1:]).sort_index()
    pddf.T.to_csv(output_file)

    print(f"Saved summarized data to {output_file}")

if __name__ == "__main__":
    main()
