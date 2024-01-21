
# assuming anaconda or miniconda is installed

# set up a python environment (pdb-tools is incuded in the environment.yml file)
# conda create --name EulerAngleEnv --file environment.yml
# conda activate EulerAngleEnv

# ---------------------------------
# Atom Selection 1
# ---------------------------------

mkdir -p atomsel1

# Principal Axes Alignment
python ../alignPrincipalAxes/alignPrincipalAxes.py --infile ./inputs/pdb2wf5.pdb \
                                                   --select "protein and resid 16:87 and name CA" \
                                                   --outfile ./atomsel1/pdb2wf5_PAaligned.pdb

# Euler Angle Calculation
python ../getEulerAngles/getEulerAngles.py --queryStructures queryStructures \
                                           --referenceFile ./atomsel1/pdb2wf5_PAaligned.pdb \
                                           --referenceChain A \
                                           --domain1 16-87 \
                                           --domain2 1-15,88-221 \
                                           --atoms CA \
                                           --outfile ./atomsel1/eulerAngles

# ---------------------------------
# Atom Selection 2
# ---------------------------------

mkdir -p atomsel2

# Principal Axes Alignment
python ../alignPrincipalAxes/alignPrincipalAxes.py --infile ./inputs/pdb2wf5.pdb \
                                                   --select "protein and resid 16:87 and name CA" \
                                                   --outfile ./atomsel2/pdb2wf5_PAaligned.pdb

# Euler Angle Calculation
python ../getEulerAngles/getEulerAngles.py --queryStructures queryStructures \
                                           --referenceFile ./atomsel2/pdb2wf5_PAaligned.pdb \
                                           --referenceChain A \
                                           --domain1 15-87 \
                                           --domain2 1-14,88-221 \
                                           --atoms CA \
                                           --outfile ./atomsel2/eulerAngles


# ---------------------------------
# Atom Selection 3
# ---------------------------------

mkdir -p atomsel3

# Principal Axes Alignment
python ../alignPrincipalAxes/alignPrincipalAxes.py --infile ./inputs/pdb2wf5.pdb \
                                                   --select "protein and resid 16:87 and name CA" \
                                                   --outfile ./atomsel3/pdb2wf5_PAaligned.pdb

# Euler Angle Calculation
python ../getEulerAngles/getEulerAngles.py --queryStructures queryStructures \
                                           --referenceFile ./atomsel3/pdb2wf5_PAaligned.pdb \
                                           --referenceChain A \
                                           --domain1 15-88 \
                                           --domain2 1-14,89-221 \
                                           --atoms CA \
                                           --outfile ./atomsel3/eulerAngles
