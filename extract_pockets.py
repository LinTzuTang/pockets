import Bio.PDB

# Standard Libraries
import os
import sys

# Biopython Libraries
from Bio.PDB import PDBParser, Select
from Bio.PDB.Residue import Residue

# RDKit Libraries
from rdkit import Chem

# Scipy Libraries
from scipy.spatial import KDTree


# xpdb.py -- extensions to Bio.PDB
# (c) 2009 Oliver Beckstein
# Relased under the same license as Biopython.
# See http://biopython.org/wiki/Reading_large_PDB_files

class SloppyStructureBuilder(Bio.PDB.StructureBuilder.StructureBuilder):
    """Cope with resSeq < 10,000 limitation by just incrementing internally.

    # Q: What's wrong here??
    #   Some atoms or residues will be missing in the data structure.
    #   WARNING: Residue (' ', 8954, ' ') redefined at line 74803.
    #   PDBConstructionException: Blank altlocs in duplicate residue SOL
    #   (' ', 8954, ' ') at line 74803.
    #
    # A: resSeq only goes to 9999 --> goes back to 0 (PDB format is not really
    #    good here)
    """

    # NOTE/TODO:
    # - H and W records are probably not handled yet (don't have examples
    #   to test)

    def __init__(self, verbose=False):
        Bio.PDB.StructureBuilder.StructureBuilder.__init__(self)
        self.max_resseq = -1
        self.verbose = verbose

    def init_residue(self, resname, field, resseq, icode):
        """Initiate a new Residue object.

        Arguments:
        o resname - string, e.g. "ASN"
        o field - hetero flag, "W" for waters, "H" for
            hetero residues, otherwise blanc.
        o resseq - int, sequence identifier
        o icode - string, insertion code

        """
        if field != " ":
            if field == "H":
                # The hetero field consists of
                # H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)

        if resseq > self.max_resseq:
            self.max_resseq = resseq

        if field == " ":
            fudged_resseq = False
            while self.chain.has_id(res_id) or resseq == 0:
                # There already is a residue with the id (field, resseq, icode)
                # resseq == 0 catches already wrapped residue numbers which
                # do not trigger the has_id() test.
                #
                # Be sloppy and just increment...
                # (This code will not leave gaps in resids... I think)
                #
                # XXX: shouldn't we also do this for hetero atoms and water??
                self.max_resseq += 1
                resseq = self.max_resseq
                res_id = (field, resseq, icode)  # use max_resseq!
                fudged_resseq = True

            if fudged_resseq and self.verbose:
                sys.stderr.write(
                    "Residues are wrapping (Residue "
                    + "('%s', %i, '%s') at line %i)."
                    % (field, resseq, icode, self.line_counter)
                    + ".... assigning new resid %d.\n" % self.max_resseq
                )
        residue = Residue(res_id, resname, self.segid)
        self.chain.add(residue)
        self.residue = residue

def read_pdb(pdb_file, name=None):
    """Read pdb file into Biopython structure. (not support residue number > 9999)"""
    if name is None:
        name = os.path.basename(pdb_file)
    parser = Bio.PDB.PDBParser(QUIET=False)
    print(f"Reading PDB file: {pdb_file}")
    return parser.get_structure(name, pdb_file)

def read_pdb_sloppy(pdb_file, name=None):
    """Read pdb file into Biopython structure."""
    if name is None:
        name = os.path.basename(pdb_file)
    sloppyparser = PDBParser(
    PERMISSIVE=True, structure_builder=SloppyStructureBuilder(),QUIET=True
    )
    return sloppyparser.get_structure(name, pdb_file)

def get_ligand_pdb(ligand_pdb):
    """Read ligand into RDKit Mol."""
    lig = Chem.MolFromPDBFile(str(ligand_pdb), removeHs=True)
    if lig is None:
        print('Failed to load ligand from', ligand_pdb)
        return None
    print(f"Ligand loaded from {ligand_pdb}")
    return lig

class PocketSelect(Select):
    """Selection class for subsetting RNA to key binding residues."""
    def __init__(self, reslist):
        self.reslist = reslist

    def accept_residue(self, residue):
        return residue in self.reslist

def get_pocket_res(rna, ligand, dist):
    """Extract residues within specified distance of ligand."""
    print("Extracting pocket residues...")
    prot_atoms = [a for a in rna.get_atoms()]
    prot_coords = [atom.get_coord() for atom in prot_atoms]

    lig_coords = []
    for i in range(ligand.GetNumAtoms()):
        pos = ligand.GetConformer().GetAtomPosition(i)
        lig_coords.append([pos.x, pos.y, pos.z])

    kd_tree = KDTree(prot_coords)
    key_pts = kd_tree.query_ball_point(lig_coords, r=dist, p=2.0)
    key_pts = set([k for l in key_pts for k in l])

    key_residues = set()
    for i in key_pts:
        atom = prot_atoms[i]
        res = atom.get_parent()
        if res.get_resname() != 'HOH':
            key_residues.add(res)
    print(f"Found {len(key_residues)} pocket residues")
    return key_residues

def write_files(rna, pocket, out_path):
    """Writes cleaned structure files for RNA and pocket."""
    io = Bio.PDB.MMCIFIO()
    io.set_structure(rna)
    io.save(out_path, PocketSelect(pocket))
    print(f"Pocket file saved to {out_path}")

def process_pdb_files(ligand_file, rna_file, out_path, dist=6.0):
    """Main processing function for a single ligand-RNA pair."""
    print(f"Processing files: Ligand - {ligand_file}, RNA - {rna_file}")
    rna = read_pdb_sloppy(rna_file)
    ligand = get_ligand_pdb(ligand_file)
    if ligand is None:
        print('Ligand processing failed for', ligand_file)
        return
    pocket_res = get_pocket_res(rna, ligand, dist)
    write_files(rna, pocket_res, out_path)

if __name__ == "__main__":
    ligand_dir = 'general_dataset/parsed_ligands_pdb'
    rna_dir_pdb = 'general_dataset/parsed_rna3db_rnas_whole_pdb'
    output_dir = 'general_dataset/pockets_6A'
    
    os.makedirs(output_dir, exist_ok=True)
    
    files = os.listdir(ligand_dir)
    total_files = len(files)
    
    for idx, ligand_file in enumerate(files, start=1):
        pdbid = ligand_file.split('_')[0]
        ligand_name = "_".join(ligand_file.split('_')[:3]).replace('.pdb', '')
        rna_file = os.path.join(rna_dir_pdb, f'{pdbid}_rna.pdb')
        
        if not os.path.exists(rna_file):
            print(f'RNA file not found for {pdbid}')
            continue
        
        output_file = os.path.join(output_dir, f'{ligand_name}_pocket.cif')
        process_pdb_files(os.path.join(ligand_dir, ligand_file), rna_file, output_file)
        
        # Print progress
        print(f"Processed {idx}/{total_files} files")
