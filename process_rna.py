import os
from pockets.rna_parser import StructureFile, ModificationHandler
from convert_cif2pdb import convert_cif_to_pdb

def process_and_convert_files(complexes_dir, output_dir_cif, output_dir_pdb):
    # Create output directories if they don't exist
    os.makedirs(output_dir_cif, exist_ok=True)
    os.makedirs(output_dir_pdb, exist_ok=True)
    
    modification_handler = ModificationHandler()

    # Process all files in the complexes_dir
    for filename in os.listdir(complexes_dir):
        if filename.endswith('.cif'):
            # Construct the full file paths
            input_path = os.path.join(complexes_dir, filename)
            output_filename_cif = f"{os.path.splitext(filename)[0]}_rna.cif"
            output_path_cif = os.path.join(output_dir_cif, output_filename_cif)
            
            # Process the structure file
            complex_structure_rna3db = StructureFile(input_path, modification_handler)
            complex_structure_rna3db.write_mmcif_chain(output_path_cif)

            # Convert CIF to PDB
            output_filename_pdb = f"{os.path.splitext(filename)[0]}_rna.pdb"
            output_path_pdb = os.path.join(output_dir_pdb, output_filename_pdb)
            convert_cif_to_pdb(output_path_cif, output_path_pdb)

            print(f"Processed and converted: {filename} -> {output_filename_cif} -> {output_filename_pdb}")

def main():
    # Define directories
    complexes_dir = 'complex_mmcif_database_non_rnasm_removed'
    output_dir_cif = 'general_dataset/parsed_rna3db_rnas_whole'
    output_dir_pdb = 'general_dataset/parsed_rna3db_rnas_whole_pdb'
    
    # Process and convert files
    process_and_convert_files(complexes_dir, output_dir_cif, output_dir_pdb)

if __name__ == "__main__":
    main()
