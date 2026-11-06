#!/usr/bin/env bash
# Example HH-suite workflow for testing homology between Influenza PA endonuclease (6FS6) and Andes L endonuclease (6Q99)
# Requirements: hh-suite (hhblits, hhmake, hhsearch) and databases (e.g., uniclust30, pdb70)
# Set these to your local DB paths
UNICLUST_DB="/path/to/uniclust30/uniclust30_2018_08/uniclust30_2018_08"
PDB70_DB="/path/to/pdb70/pdb70"

mkdir -p hhsuite
cd hhsuite || exit 1

# Inputs from this project
cp ../analysis/sequences/6FS6_chainE.fasta .
cp ../analysis/sequences/6Q99_chainA.fasta .

# Build MSAs and HHMs
hhblits -i 6FS6_chainE.fasta -oa3m 6FS6.a3m -n 3 -d "$UNICLUST_DB" -cpu 4
hhmake -i 6FS6.a3m -o 6FS6.hhm

hhblits -i 6Q99_chainA.fasta -oa3m 6Q99.a3m -n 3 -d "$UNICLUST_DB" -cpu 4
hhmake -i 6Q99.a3m -o 6Q99.hhm

# Direct profile–profile comparison
hhsearch -i 6FS6.hhm -d 6Q99.hhm -o FS6_vs_Q99.hhr -cpu 4
hhsearch -i 6Q99.hhm -d 6FS6.hhm -o Q99_vs_FS6.hhr -cpu 4

# Optional: search against PDB70 for context
hhsearch -i 6FS6.hhm -d "$PDB70_DB" -o FS6_vs_pdb70.hhr -cpu 4
hhsearch -i 6Q99.hhm -d "$PDB70_DB" -o Q99_vs_pdb70.hhr -cpu 4

echo "Done. Inspect *.hhr for probability and E-values. Significant matches support homology (divergent evolution)."
