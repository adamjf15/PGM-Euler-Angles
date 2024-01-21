import MDAnalysis as mda
import argparse
import warnings

warnings.filterwarnings(
    action="ignore",
    module="MDAnalysis",
    message="Found no information for attr: 'formalcharges'"
)

# --------------------
# Functions & Utils
# --------------------

def align_to_principal_axes(ag):

    ag.universe.atoms.translate(-ag.center_of_mass())
    ag.universe.atoms.rotate(ag.principal_axes())

    return None

# --------------------
# Argument Parsing
# --------------------

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile", required=True, type=str)
parser.add_argument("-s", "--select", required=False, default="protein", type=str)
parser.add_argument("-o", "--outfile", required=True, type=str)
args = parser.parse_args()

# --------------------
# Main
# --------------------

u = mda.Universe(args.infile)

ag = u.select_atoms(args.select)

align_to_principal_axes(ag)

with mda.Writer(args.outfile) as f:
    f.write(u.atoms)