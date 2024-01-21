import warnings
import numpy as np
import pandas as pd
import argparse
from scipy.spatial.transform import Rotation as R
from pdbfixer.pdbfixer import substitutions
from Bio import Align
from Bio import BiopythonDeprecationWarning

# ---------------------------
# Warning Filters
# ---------------------------

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=BiopythonDeprecationWarning)
    import MDAnalysis as mda
    from MDAnalysis.core.groups import AtomGroup, ResidueGroup
    from MDAnalysis.analysis import align

warnings.filterwarnings(
    action="ignore",
    module="MDAnalysis",
    message="Failed to guess the mass for the following atom types"
)

# ---------------------------
# Catalogued Residues
# ---------------------------

# standard residues
standard_residues = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                     'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                     "PHD": 'D', "FTR": "W"
                    }

# known non-standard residues from OpenMM pdb-fixer module
nonstandard_residues = {key: standard_residues[substitutions[key]] for key in substitutions.keys()}

# updates in-place
known_residues = {**standard_residues, **nonstandard_residues} 

# ---------------------------
# Functions / Utils
# ---------------------------

def euler_angles(QueryFile, QueryChain, ReferenceFile, ReferenceChain, Domain1Resids, Domain2Resids, atoms):

    # initialise MDA objects
    Query = mda.Universe(QueryFile)
    Reference = mda.Universe(ReferenceFile)

    # create residue selections
    QueryRG = Query.select_atoms(f'protein and chainID {QueryChain}').residues
    ReferenceRG = Reference.select_atoms(f'protein and chainID {ReferenceChain}').residues

    # remove unrecognised residues
    QueryRG = remove_unrecognised_residues(QueryRG)
    ReferenceRG = remove_unrecognised_residues(ReferenceRG)

    # Get matched residue selections via sequence alignment
    QueryRG, ReferenceRG = match_residues(QueryRG, ReferenceRG)

    # Domain1 Residue Selection
    idx1 = [i for i, r in enumerate(ReferenceRG) if r.resid in Domain1Resids]
    ReferenceRG1 = ReferenceRG[idx1]
    QueryRG1 = QueryRG[idx1]
    
    # Domain2 Residue Selection
    idx2 = [i for i, r in enumerate(ReferenceRG) if r.resid in Domain2Resids]
    ReferenceRG2 = ReferenceRG[idx2]
    QueryRG2 = QueryRG[idx2]

    # Get matched atom selections from matched residues
    QueryAG1, ReferenceAG1 = match_atoms_from_matched_residue_groups(QueryRG1, ReferenceRG1, atoms)
    QueryAG2, ReferenceAG2 = match_atoms_from_matched_residue_groups(QueryRG2, ReferenceRG2, atoms)

    # Align to Domain 2
    align.alignto(QueryAG2, ReferenceAG2)

    # Alignment to Domain 1
    rot, rmsd = align.rotation_matrix(QueryAG1.positions, ReferenceAG1.positions)

    # calculate euler angles from second alignment
    yaw, pitch, roll = R.from_matrix(rot.T).as_euler('XYZ', degrees=True)

    # calculate screw-axis angle from second alignment
    theta = np.degrees(np.arccos((np.trace(rot) - 1) / 2))

    # return angles
    return np.array([yaw, pitch, roll, theta])


def remove_unrecognised_residues(rg):

    return ResidueGroup([r for r in rg if r.resname in known_residues.keys()])


def match_residues(rg1, rg2):

    # get protein sequences
    pepseq1 = get_peptide_sequece(rg1)
    pepseq2 = get_peptide_sequece(rg2)

    # align protein sequences
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(pepseq1, pepseq2)
    alignment = alignments[0]

    # get indicies of equivalent residues
    mask = np.all(alignment.indices != -1, axis=0)
    inds1, inds2 = alignment.indices[:, mask]

    # get equivalent residues
    rg1 = rg1[inds1]
    rg2 = rg2[inds2]

    # return equivalent residues
    return rg1, rg2


def get_peptide_sequece(rg):

    return "".join(known_residues[x] for x in rg.resnames)


def match_atoms_from_matched_residue_groups(rg1, rg2, atoms):

    # initialise empty atomgroups
    ag1 = AtomGroup([], rg1.universe)
    ag2 = AtomGroup([], rg2.universe)

    for r1, r2 in zip(rg1, rg2):

        r1_ag, r2_ag = match_atoms_from_matched_residues(r1, r2, atoms=atoms)
        ag1 = ag1.concatenate(r1_ag)
        ag2 = ag2.concatenate(r2_ag)

    return ag1, ag2


def match_atoms_from_matched_residues(r1, r2, atoms):

    r1_atoms = [a for a in r1.atoms if (a.name in atoms) and (a.name in r2.atoms.names)]
    r1_atoms = sorted(r1_atoms, key=lambda a: a.name)
    r1_atoms = AtomGroup(r1_atoms)
    r1_atoms = remove_duplicate_atoms(r1_atoms)
    
    r2_atoms = [a for a in r2.atoms if (a.name in atoms) and (a.name in r1.atoms.names)]
    r2_atoms = sorted(r2_atoms, key=lambda a: a.name)
    r2_atoms = AtomGroup(r2_atoms)
    r2_atoms = remove_duplicate_atoms(r2_atoms)

    assert (len(r1_atoms) == len(r2_atoms))
    assert all(r1_atoms.names == r2_atoms.names)

    return r1_atoms, r2_atoms


def remove_duplicate_atoms(ag):

    ids = [f'{a.resname}-{a.resid}-{a.name}' for a in ag]
    _, idx = np.unique(ids, return_index=True)

    return ag[idx]


# ---------------------------
# Batch Calculations
# ---------------------------


def batch_euler_angles(QueryFiles, QueryChains, ReferenceFile, ReferenceChain, Domain1Resids, Domain2Resids, atoms):

    data = []
    for QueryFile, QueryChain in zip(QueryFiles, QueryChains):
        StructureName = f'{QueryFile.split("/")[-1].split(".")[0]}_{QueryChain}'
        yaw, pitch, roll, theta = euler_angles(QueryFile, QueryChain, ReferenceFile, ReferenceChain, Domain1Resids, Domain2Resids, atoms)
        data.append({"structure": StructureName,
                     "yaw": yaw,
                     "pitch": -pitch,
                     "roll": roll,
                     "theta": theta})

        output_string = f'{StructureName:<20}yaw: {yaw:5.1f}  pitch: {-pitch:5.1f}  roll: {roll:5.1f}  theta: {theta:5.1f}'
        print(output_string)
        

    df = pd.DataFrame(data)

    return df

# ---------------------------
# Wrting Output
# ---------------------------

def write_batch_euler_angles(outfile, EulerAngles):

    f = open(outfile, "w")
    f.write(EulerAngles.to_string(index=None, float_format="%8.3f"))
    f.close()

# ---------------------------
# Parsing Functions
# ---------------------------


def parse_range_string(range_string):
    
    result = []
    
    ranges = range_string.split(',')
    
    for r in ranges:
        if '-' in r:
            start, end = map(int, r.split('-'))
            result.extend(range(start, end + 1))
        else:
            result.append(int(r))
    
    return result

def parse_input_file(QueryInputs):

    with open(QueryInputs, "r") as f:
        contents = f.read()

    lines = contents.split("\n")
    
    QueryFiles = [line.split()[0] for line in lines]
    QueryChains = [line.split()[1] for line in lines]

    return QueryFiles, QueryChains

# ---------------------------
# Config
# ---------------------------

VERSION = "v0.1.8"

# ---------------------------
# Argument Parsing
# ---------------------------

program_description="""
A python program for calculating Euler angles quantifying the
orientation of one domain relative to another in a protein
structure.
"""

# create argument parser
parser = argparse.ArgumentParser(program_description)

# query structures
parser.add_argument("--queryStructures", required=True, type=str)

# reference structure
parser.add_argument("--referenceFile", required=True, type=str)
parser.add_argument("--referenceChain", required=True, type=str)
parser.add_argument("--domain1", required=True, type=str)
parser.add_argument("--domain2", required=True, type=str)
parser.add_argument("--atoms", required=True, type=str, nargs="*")

# ouput
parser.add_argument("-o", "--outfile", required=True, type=str)

# parsing
args = parser.parse_args()


# ---------------------------
# Main
# ---------------------------

# parse query structures
QueryFiles, QueryChains = parse_input_file(args.queryStructures)

# run alignments & calculations
EulerAngles = batch_euler_angles(
    QueryFiles=QueryFiles, 
    QueryChains=QueryChains,
    ReferenceFile=args.referenceFile,
    ReferenceChain=args.referenceChain,
    Domain1Resids=parse_range_string(args.domain1),
    Domain2Resids=parse_range_string(args.domain2),
    atoms=args.atoms,
)

# write output
write_batch_euler_angles(
    outfile=args.outfile,
    EulerAngles=EulerAngles,
)

