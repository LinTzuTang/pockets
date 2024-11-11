import os
import gemmi
import pandas as pd
import argparse


# Function to convert a single CIF to PDB using gemmi
def convert_cif_to_pdb(cif_file, pdb_file):
    doc = gemmi.cif.read_file(cif_file)
    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)
    changes = []

    # Check and modify chain names and residue sequence numbers if they are too long
    for model in structure:
        for chain in model:
            original_name = chain.name
            new_chain_name = original_name
            if len(chain.name) > 1:
                new_chain_name = chain.name[0]  # Shorten chain name to first character
                chain.name = new_chain_name
                print(
                    f"Chain name {original_name} in {cif_file} is too long, shortened to {chain.name}."
                )
                changes.append(
                    {
                        "cif_file": os.path.basename(cif_file),
                        "chain_name": original_name,
                        "new_chain_name": new_chain_name,
                    }
                )
            for residue in chain:
                original_seq_id = residue.seqid.num
                new_seq_id = original_seq_id
                if residue.seqid.num > 9999:
                    new_seq_id = 9999  # Set residue sequence number to 9999
                    residue.seqid.num = new_seq_id
                    print(f"Residue sequence number {original_seq_id} in {cif_file} is too large, set to {residue.seqid.num}.")

            # for residue in chain:
            #     if residue.seqid.num > 9999:
            #         original_seq_id = residue.seqid.num
            #         new_seq_id = hex(original_seq_id)  # Convert to hexadecimal
            #         residue.seqid.num = int(new_seq_id, 16)
            #         print(f"Residue sequence number {original_seq_id} in {cif_file} is too large, converted to {residue.seqid.num}.")

    structure.write_pdb(pdb_file)
    return changes


# Function to convert all CIF files in a directory
def convert_all_cif_in_dir(input_cif_dir, output_pdb_dir):
    # Ensure the output directory exists
    os.makedirs(output_pdb_dir, exist_ok=True)

    # DataFrame to store the changes
    changes_df = pd.DataFrame(columns=["cif_file", "chain_name", "new_chain_name"])

    # Loop through all CIF files in the directory
    for cif_filename in os.listdir(input_cif_dir):
        if cif_filename.endswith(".cif"):
            cif_file = os.path.join(input_cif_dir, cif_filename)
            pdb_filename = cif_filename.replace(".cif", ".pdb")
            pdb_file = os.path.join(output_pdb_dir, pdb_filename)

            # Convert CIF to PDB and record changes
            changes = convert_cif_to_pdb(cif_file, pdb_file)
            changes_df = changes_df.append(changes, ignore_index=True)
            print(f"Converted {cif_file} to {pdb_file}")

    # Save the changes to a CSV file in the output directory
    changes_csv_path = os.path.join(output_pdb_dir, "cif_to_pdb_changes.csv")
    changes_df.to_csv(changes_csv_path, index=False)


def main():
    parser = argparse.ArgumentParser(description="Convert CIF files to PDB format.")
    parser.add_argument(
        "input_cif_dir", type=str, help="Directory containing CIF files"
    )
    parser.add_argument(
        "output_pdb_dir", type=str, help="Directory to save the converted PDB files"
    )

    args = parser.parse_args()

    convert_all_cif_in_dir(args.input_cif_dir, args.output_pdb_dir)


if __name__ == "__main__":
    # Uncomment one of the following two blocks to use the respective mode:

    # Mode 1: Use command-line arguments
    # main()

    # Mode 2: Use hardcoded directories
    cif_dir = "../hariboss_process/hariboss_20240621_error_revised_result/parsed_rna3db_rnas"
    output_dir = "../hariboss_process/hariboss_20240621_error_revised_result/parsed_rna3db_rnas_pdb"
    convert_all_cif_in_dir(cif_dir, output_dir)
