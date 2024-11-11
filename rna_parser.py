import dataclasses
from Bio import PDB

import os
import json
from typing import Sequence
from collections import defaultdict


class ModificationHandler:
    JSON_PATH_DEFAULTS = [
        "modifications_cache.json",
        "data/modifications_cache.json",
    ]

    def __init__(self, json_path: str = None):
        """Used for converting `three_letter_code`s to `one_letter_code`s, including modifications.

        Args:
            json_path (str, optional): Path to `modifications_cache.json`, generated from a Chemical Component
                Dictionary .cif file. If not provided, `JSON_PATH_DEFAULTS` are checked for a valid path.

        Attributes:
            JSON_PATH_DEFAULTS (Sequence[str]): List of paths checked for a `modifications_cache.json` file when
                `json_path` is not initialised.
        """
        # try to see we can find the modifications_cache.json in any of the default paths
        if json_path is None:
            for path in self.JSON_PATH_DEFAULTS:
                if os.path.isfile(path):
                    json_path = path
                    break

        # if we still haven't found it, we need to throw an exception
        if json_path is None:
            raise FileNotFoundError(
                "Could not find `modifications_cache.json` in JSON_PATH_DEFAULTS."
            )

        with open(json_path) as f:
            self.modifications = json.load(f)

    def is_rna(self, three_letter_code: str) -> bool:
        """Check if `three_letter_code` is RNA nucleic acid.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
            bool: True if `three_letter_code` is typically an RNA nucleic acid, False otherwise.

        """
        return three_letter_code in self.modifications["rna"]

    def is_protein(self, three_letter_code: str) -> bool:
        """Check if `three_letter_code` is an amino acid.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
            bool: True if `three_letter_code` is an amino acid.

        """
        return three_letter_code in self.modifications["protein"]

    def protein_letters_3to1(self, three_letter_code: str) -> str:
        """Convert amino acid `three_letter_code` to `one_letter_code`.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
           str: one_letter_code of an amino acid, "X" if cannot be found.

        """
        return self.modifications["protein"].get(three_letter_code, "X")

    def rna_letters_3to1(self, three_letter_code: str) -> str:
        """Convert RNA nucleic acid `three_letter_code` to `one_letter_code`.

        Args:
            three_letter_code (str): Three letter code to check.

        Returns:
           str: one_letter_code of RNA nucleic acid, "N" if cannot be found.
        """
        return self.modifications["rna"].get(three_letter_code, "N")


class Residue:
    """Data class wrapping individual residues."""

    def __init__(
        self,
        three_letter_code: str,
        one_letter_code: str,
        index: int,
    ):
        self.three_letter_code = three_letter_code
        self.one_letter_code = one_letter_code
        self.index = index
        self.atoms = {}

    @property
    def code(self) -> str:
        return self.one_letter_code

    @property
    def is_missing(self) -> bool:
        return not len(self.atoms) > 0

    def __repr__(self):
        return (
            f"Residue(code={self.code}, three_letter_code={self.three_letter_code}, "
            f"index={self.index}, is_missing={self.is_missing})"
        )


class Chain:
    """Data class wrapping chains. Contains a list of Residues."""

    def __init__(self, author_id: str = None):
        self.author_id = author_id
        self.residues = []

    def __iter__(self):
        if len(self) == 0:
            return None
        return iter(self.residues)

    def __getitem__(self, idx):
        if len(self) == 0:
            return None
        return self.residues[idx]

    def __len__(self):
        return len(self.residues)

    @property
    def has_atoms(self):
        return any([not res.is_missing for res in self])

    def add_residue(self, res: Residue):
        """Add a residue to the chain.

        NOTE: Residues indices should be added in increasing order.
              If two conecutive residues are added with the same index, we
              ignore the second one. See 4x4t_G for the motivation.
              We also infer missing residues as `N`s if there is a gap in the
              insertion.

        Args:
            res (Residue): Residue to add to the chain.
        """
        if len(self) == 0 or self.residues[-1].index == res.index - 1:
            self.residues.append(res)
        elif self.residues[-1].index < res.index - 1:
            for idx in range(self.residues[-1].index + 1, res.index):
                self.residues.append(Residue("N", "N", idx))
            self.residues.append(res)
        elif self.residues[-1].index == res.index:
            pass
        else:
            raise ValueError(f"Cannot add residues out of order.")

    @property
    def sequence(self):
        return "".join(i.code for i in self.residues)

    def __repr__(self):
        return f"Chain(author_id={self.author_id}, len={len(self)})"

    def __str__(self):
        max_line_length = 120
        idx_steps = 50

        numbers_str = [" " for i in self.sequence]
        lbl = list(str(len(self.sequence)))
        numbers_str[-len(lbl) :] = lbl
        numbers_str[0] = "1"

        for idx in range(idx_steps, len(self.sequence), idx_steps):
            lbl = list(str(idx + 1))
            numbers_str[idx - len(lbl) : idx] = lbl

        numbers_str = "".join(numbers_str)

        S = ""
        for idx in range(0, len(self.sequence), max_line_length):
            S += numbers_str[idx : idx + max_line_length] + "\n"
            S += self.sequence[idx : idx + max_line_length] + "\n"
            S += "\n"

        return S


class mmCIFParser:
    def __init__(
        self,
        path: str,
        modification_handler: ModificationHandler,
    ):
        self.path = path
        self.molecule_type = "RNA"
        self.nmr_resolution = None
        self.include_atoms = True

        self.letters_3to1 = lambda x: modification_handler.rna_letters_3to1(x)

        self.parsed_info = PDB.MMCIF2Dict.MMCIF2Dict(self.path)

    @property
    def pdb_id(self):
        return self.parsed_info["_entry.id"][0].lower()

    @property
    def release_date(self):
        # Prefer to get the date from the earliest revision date
        if "_pdbx_audit_revision_history.revision_date" in self.parsed_info:
            return min(self.parsed_info["_pdbx_audit_revision_history.revision_date"])
        # Use deposition date if there are no revisions
        return self.parsed_info["_pdbx_database_status.recvd_initial_deposition_date"]

    @property
    def resolution(self):
        resolutions = []
        for res_key in [
            "_refine.ls_d_res_high",
            "_em_3d_reconstruction.resolution",
            "_reflns.d_resolution_high",
        ]:
            if res_key in self.parsed_info:
                if self.parsed_info[res_key][0] not in ".?":
                    resolutions.append(float(self.parsed_info[res_key][0]))

        # If we have an NMR structure and we overwrite default NMR resolution
        if self.structure_method == "solution nmr" and self.nmr_resolution is not None:
            return self.nmr_resolution

        if len(resolutions) == 0:
            return float("inf")

        return max(resolutions)

    @property
    def structure_method(self):
        return ",".join(self.parsed_info["_exptl.method"]).lower()

    @dataclasses.dataclass
    class _AtomSite:
        atom_id: str
        three_letter_code: str
        author_chain_id: str
        entity_id: str
        author_seq_num: int
        mmcif_seq_num: int
        insertion_code: str
        hetatm_atom: str
        alt_id: str
        x: float
        y: float
        z: float
        model_num: int

    @staticmethod
    def _get_atom_sites(parsed_info: PDB.MMCIF2Dict):
        return [
            mmCIFParser._AtomSite(*site)
            for site in zip(
                parsed_info["_atom_site.label_atom_id"],  # atom name
                parsed_info["_atom_site.label_comp_id"],  # residue name
                parsed_info["_atom_site.auth_asym_id"],  # author chain
                parsed_info["_atom_site.label_entity_id"],  # entity id
                parsed_info["_atom_site.auth_seq_id"],  # author_seq num
                parsed_info["_atom_site.label_seq_id"],  # mmcif_seq_num
                parsed_info["_atom_site.pdbx_PDB_ins_code"],  # insertion code?
                parsed_info["_atom_site.group_PDB"],  # hetatm_atom
                parsed_info[
                    "_atom_site.label_alt_id"
                ],  # alternative conformation identifier
                parsed_info["_atom_site.Cartn_x"],  # x
                parsed_info["_atom_site.Cartn_y"],  # y
                parsed_info["_atom_site.Cartn_z"],  # z
                parsed_info["_atom_site.pdbx_PDB_model_num"],  # model
            )
        ]

    @property
    def chains(self):
        # No SEQRES chains in this file
        if "_entity_poly_seq.entity_id" not in self.parsed_info:
            return {}

        # Get the mapping from entity_id to internal mmcif_chain_id
        mmcif_chain_to_entity_id = {
            mmcif_chain_id: entity_id
            for mmcif_chain_id, entity_id in zip(
                self.parsed_info["_struct_asym.id"],
                self.parsed_info["_struct_asym.entity_id"],
            )
        }

        # Create mapping that maps entity_id to author_chain_id
        id_map = defaultdict(set)
        for author_chain_id, mmcif_chain_id in zip(
            self.parsed_info["_atom_site.auth_asym_id"],
            self.parsed_info["_atom_site.label_asym_id"],
        ):
            k = mmcif_chain_to_entity_id[mmcif_chain_id]
            id_map[k].add(author_chain_id)

        # Get the chem_comp type for each mon_id
        chem_comp_type = {
            mon_id: comp_type
            for mon_id, comp_type in zip(
                self.parsed_info["_chem_comp.id"], self.parsed_info["_chem_comp.type"]
            )
        }

        # Parse full chains from "seqres"
        chains_full = defaultdict(Chain)
        for entity_id, mon_id, idx in zip(
            self.parsed_info["_entity_poly_seq.entity_id"],
            self.parsed_info["_entity_poly_seq.mon_id"],
            self.parsed_info["_entity_poly_seq.num"],
        ):
            for author_id in id_map[entity_id]:
                chains_full[author_id].author_id = author_id
                chains_full[author_id].add_residue(
                    Residue(
                        three_letter_code=mon_id,
                        one_letter_code=self.letters_3to1(mon_id),
                        index=int(idx) - 1,
                    )
                )

        # Keep only chains that contain at least one 'self.molecule_type'
        chains = {}
        for author_chain_id, chain_data in chains_full.items():
            if any(
                [
                    self.molecule_type in chem_comp_type[i.three_letter_code]
                    for i in chain_data.residues
                ]
            ):
                chains[author_chain_id] = chain_data

        # Find starting index of relevant chains
        seq_start_num = {
            author_chain_id: min([res.index for res in chain])
            for author_chain_id, chain in chains.items()
        }

        # Keep track of alt conformations so we only take one
        chain_alt_id = {}

        # Iterate through atom sites to get coordinates if required
        if self.include_atoms:
            for site in mmCIFParser._get_atom_sites(self.parsed_info):
                # Get RNA only

                if site.model_num != "1":  # Get model 1 only
                    continue
                # If not a relevant chain, don't care about atoms
                if site.author_chain_id not in chains:
                    continue

                # Handle missing sequence numbers
                if site.mmcif_seq_num == ".":
                    continue

                # Make sure we always take the same alt conformation
                if site.author_chain_id not in chain_alt_id:
                    chain_alt_id[site.author_chain_id] = site.alt_id
                if site.alt_id != chain_alt_id[site.author_chain_id]:
                    continue

                # The idx is just the mmcif_seq_num - seq_start_num
                seq_idx = (
                    int(site.mmcif_seq_num) - seq_start_num[site.author_chain_id] - 1
                )

                # Handle potential mismatches
                if (
                    site.three_letter_code
                    != chains[site.author_chain_id][seq_idx].three_letter_code
                ):
                    print(
                        f"WARNING: found a mismatch in {self.pdb_id}_{site.author_chain_id} ({site.entity_id}) at position {seq_idx}, "
                        f"(entity_poly_seq: {chains[site.author_chain_id][seq_idx].three_letter_code} atom_site: {site.three_letter_code})"
                    )
                    continue

                # Add atom coordinates
                chains[site.author_chain_id][seq_idx].atoms[site.atom_id] = tuple(
                    map(float, (site.x, site.y, site.z))
                )

        return chains


class StructureFile:
    def __init__(
        self,
        path: str,
        modification_handler: ModificationHandler,
    ):
        # determine which parser to use
        if path.lower().endswith((".cif", ".mmcif")):
            file_parser = mmCIFParser
        elif path.lower().endswith(".pdb"):
            raise NotImplementedError(
                "Unable to parse PDB files. Please use PDBx/mmCIF."
            )
        else:
            raise ValueError(
                f"The extension {path.split('.')[-1].lower()} is not supported."
            )

        # make the parser with default molecule_type as "RNA" and include_atoms=True
        parser = file_parser(path, modification_handler)

        # use parser to get attributes
        self.pdb_id = parser.pdb_id
        self.release_date = parser.release_date
        self.resolution = parser.resolution
        self.structure_method = parser.structure_method
        self.chains = parser.chains
        self.path = parser.path

    def __repr__(self):
        return (
            f"StructureFile(pdb_id={self.pdb_id}, chains={self.chains.keys()}, "
            f"resolution={self.resolution}, release_date={self.release_date}, structure_method={self.structure_method})"
        )

    def __getitem__(self, idx):
        return self.chains[idx]

    def __iter__(self):
        return iter(self.chains.values())

    def get_path(self):
        return self.path

    @staticmethod
    def _gen_mmcif_loop_str(name: str, headers: Sequence[str], values: Sequence[tuple]):
        s = "#\nloop_\n"
        for header in headers:
            s += f"_{name}.{header}\n"

        max_widths = {k: 0 for k in headers}
        for V in values:
            for k, v in zip(headers, V):
                max_widths[k] = max(max_widths[k], len(str(v)))

        for V in values:
            row = ""
            for k, v in zip(headers, V):
                row += f"{str(v):<{max_widths[k]}} "
            s += row + "\n"

        return s

    def write_mmcif_chain(self, output_path):
        # Initialize containers for all chains
        entity_poly_seq_data = []
        atom_site_data = []

        # Iterate over all chains
        for chain_id, chain in self.chains.items():
            if not chain.has_atoms:
                raise ValueError(
                    f"Did not find any atoms for chain {chain_id}. Did you set `include_atoms=True`?"
                )

            # Extract needed info for the current chain
            for i, res in enumerate(chain):
                entity_poly_seq_data.append((1, res.index + 1, res.code, "n"))
                for idx, (atom_name, atom_coords) in enumerate(res.atoms.items()):
                    x, y, z = atom_coords
                    atom_site_data.append(
                        (
                            "ATOM",
                            idx + 1,
                            atom_name[0],
                            atom_name,
                            ".",
                            res.code,
                            chain_id,
                            "?",
                            i + 1,
                            "?",
                            x,
                            y,
                            z,
                            1.0,
                            0.0,
                            "?",
                            i + 1,
                            res.code,
                            chain_id,
                            atom_name,
                            1,
                        )
                    )

        # Build required strings
        header_str = (
            f"# generated by rna3db\n"
            f"#\n"
            f"data_{self.pdb_id}\n"
            f"_entry.id {self.pdb_id}\n"
            f"_pdbx_database_status.recvd_initial_deposition_date {self.release_date}\n"
            f"_exptl.method '{self.structure_method.upper()}'\n"
            f"_reflns.d_resolution_high {self.resolution}\n"
        )

        struct_asym_str = StructureFile._gen_mmcif_loop_str(
            "_struct_asym",
            [
                "id",
                "pdbx_blank_PDB_chainid_flag",
                "pdbx_modified",
                "entity_id",
                "details",
            ],
            [(chain_id, "N", "N", 1, "?") for chain_id in self.chains],
        )

        chem_comp_str = StructureFile._gen_mmcif_loop_str(
            "_chem_comp",
            [
                "id",
                "type",
                "mon_nstd_flag",
                "pdbx_synonyms",
                "formula",
                "formula_weight",
            ],
            [
                (
                    "A",
                    "'RNA linking'",
                    "y",
                    '"ADENOSINE-5\'-MONOPHOSPHATE"',
                    "?",
                    "'C10 H14 N5 O7 P'",
                    347.221,
                ),
                (
                    "C",
                    "'RNA linking'",
                    "y",
                    '"CYTIDINE-5\'-MONOPHOSPHATE"',
                    "?",
                    "'C9 H14 N3 O8 P'",
                    323.197,
                ),
                (
                    "G",
                    "'RNA linking'",
                    "y",
                    '"GUANOSINE-5\'-MONOPHOSPHATE"',
                    "?",
                    "'C9 H13 N2 O9 P'",
                    363.221,
                ),
                (
                    "U",
                    "'RNA linking'",
                    "y",
                    '"URIDINE-5\'-MONOPHOSPHATE"',
                    "?",
                    "'C9 H13 N2 O9 P'",
                    324.181,
                ),
                ("T", "'RNA linking'", "y", '"T"', "?", "''", 0),
                ("N", "'RNA linking'", "y", '"N"', "?", "''", 0),
            ],
        )

        entity_poly_seq_str = StructureFile._gen_mmcif_loop_str(
            "entity_poly_seq",
            [
                "entity_id",
                "num",
                "mon_id",
                "heter",
            ],
            entity_poly_seq_data,
        )

        atom_site_str = StructureFile._gen_mmcif_loop_str(
            "atom_site",
            [
                "group_PDB",
                "id",
                "type_symbol",
                "label_atom_id",
                "label_alt_id",
                "label_comp_id",
                "label_asym_id",
                "label_entity_id",
                "label_seq_id",
                "pdbx_PDB_ins_code",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "occupancy",
                "B_iso_or_equiv",
                "pdbx_formal_charge",
                "auth_seq_id",
                "auth_comp_id",
                "auth_asym_id",
                "auth_atom_id",
                "pdbx_PDB_model_num",
            ],
            atom_site_data,
        )

        # Write to file
        with open(output_path, "w") as f:
            f.write(header_str)
            f.write(struct_asym_str)
            f.write(chem_comp_str)
            f.write(entity_poly_seq_str)
            f.write(atom_site_str)