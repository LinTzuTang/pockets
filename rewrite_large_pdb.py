import os
import sys
import argparse
import Bio.PDB
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Residue import Residue


class SloppyStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    """Handles the resSeq < 10,000 limitation by incrementing internally."""

    def __init__(self, verbose=False):
        super().__init__()
        self.max_resseq = -1
        self.verbose = verbose

    def init_residue(self, resname, field, resseq, icode):
        if field != " ":
            field = "H_" + resname if field == "H" else field
        res_id = (field, resseq, icode)
        if resseq > self.max_resseq:
            self.max_resseq = resseq

        if field == " ":
            fudged_resseq = False
            while self.chain.has_id(res_id) or resseq == 0:
                self.max_resseq += 1
                resseq = self.max_resseq
                res_id = (field, resseq, icode)
                fudged_resseq = True

            if fudged_resseq and self.verbose:
                sys.stderr.write(
                    f"Residues are wrapping (Residue ('{field}', {resseq}, '{icode}') "
                    f"at line {self.line_counter})... assigning new resid {self.max_resseq}.\n"
                )
        residue = Residue(res_id, resname, self.segid)
        self.chain.add(residue)
        self.residue = residue


class SloppyPDBIO(PDBIO):
    """PDBIO class for large pdb files with atom and resSeq wrapping."""

    _ATOM_FORMAT_STRING = (
        "%s%5i %-4s%c%3s %c%4i%c   " + "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"
    )

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname, resseq, icode, chain_id, element="  ", charge="  "):
        record_type = "HETATM" if hetfield != " " else "ATOM  "
        name = atom.get_fullname()
        altloc = atom.get_altloc()
        x, y, z = atom.get_coord()
        bfactor = atom.get_bfactor()
        occupancy = atom.get_occupancy()
        args = (
            record_type,
            atom_number % 100000,
            name,
            altloc,
            resname,
            chain_id,
            resseq % 10000,
            icode,
            x,
            y,
            z,
            occupancy,
            bfactor,
            segid,
            element,
            charge,
        )
        return self._ATOM_FORMAT_STRING % args


def rewrite_big_pdb(pdb_file, output_dir, name=None):
    """Read PDB file and save to output directory."""
    if name is None:
        name = os.path.basename(pdb_file)
    parser = PDBParser(PERMISSIVE=True, structure_builder=SloppyStructureBuilder())
    structure = parser.get_structure(name, pdb_file)

    # Write to output directory
    io = SloppyPDBIO()
    io.set_structure(structure)
    io.save(os.path.join(output_dir, name))


def process_all_pdbs(input_dir, output_dir):
    """Process all PDB files in input_dir and save to output_dir."""
    os.makedirs(output_dir, exist_ok=True)
    for pdb_file in os.listdir(input_dir):
        if pdb_file.endswith('.pdb'):
            pdb_path = os.path.join(input_dir, pdb_file)
            rewrite_big_pdb(pdb_path, output_dir)
            print(f"Processed and saved {pdb_file} to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Process large PDB files in a directory.")
    parser.add_argument("input_dir", help="Directory containing PDB files to process")
    parser.add_argument("output_dir", help="Directory to save processed PDB files")

    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        sys.exit(1)

    process_all_pdbs(input_dir, output_dir)


if __name__ == "__main__":
    main()
