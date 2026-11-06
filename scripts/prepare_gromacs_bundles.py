import os
from pathlib import Path
import shutil
import pandas as pd

BASE = Path('.')
DOCK_DIR = BASE / 'docking' / 'ANDV_6Q99_pocket1_metal_modified'
MDROOT = BASE / 'md_runs'
OUTROOT = BASE / 'gmx_prep'
TOP_N = 10

REDOCK_DETAILED = BASE / 'potential drugs' / 'redock_metal_modified_scores_detailed.csv'

minim_mdp = """integrator  = steep
nsteps      = 50000
emtol       = 500.0
emstep      = 0.01
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
constraints = none
"""

nvt_mdp = """integrator    = md
dt            = 0.002
nsteps        = 50000   ; 100 ps
nstxout       = 500
nstvout       = 500
nstenergy     = 250
nstlog        = 250
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
vdwtype       = Cut-off
tcoupl        = V-rescale
tc-grps       = Protein_LIG Water_and_ions
tau_t         = 0.5       0.5
ref_t         = 300       300
pcoupl        = no
constraints   = h-bonds
constraint_algorithm = lincs
lincs_iter    = 1
lincs_order   = 4
dispCorr      = EnerPres
pbc           = xyz
"""

npt_mdp = """integrator    = md
dt            = 0.002
nsteps        = 500000  ; 1 ns
nstxout       = 1000
nstvout       = 1000
nstenergy     = 500
nstlog        = 500
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
vdwtype       = Cut-off
tcoupl        = V-rescale
tc-grps       = Protein_LIG Water_and_ions
tau_t         = 0.5       0.5
ref_t         = 300       300
pcoupl        = Parrinello-Rahman
pcoupltype    = isotropic
tau_p         = 5.0
ref_p         = 1.0
compressibility = 4.5e-5
constraints   = h-bonds
constraint_algorithm = lincs
lincs_iter    = 1
lincs_order   = 4
dispCorr      = EnerPres
pbc           = xyz
"""

md_mdp = """integrator    = md
dt            = 0.002
nsteps        = 5000000 ; 10 ns
nstxout       = 2500
nstvout       = 2500
nstenergy     = 1000
nstlog        = 1000
nstxout-compressed = 500
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
vdwtype       = Cut-off
tcoupl        = V-rescale
tc-grps       = Protein_LIG Water_and_ions
tau_t         = 0.5       0.5
ref_t         = 300       300
pcoupl        = Parrinello-Rahman
pcoupltype    = isotropic
tau_p         = 5.0
ref_p         = 1.0
compressibility = 4.5e-5
constraints   = h-bonds
constraint_algorithm = lincs
lincs_iter    = 1
lincs_order   = 4
dispCorr      = EnerPres
pbc           = xyz
"""

readme_tpl = """GROMACS local run instructions (per-complex)

This bundle contains:
 - complex.pdb           : Protein (prepared, includes Mn2+) + docked ligand pose (MOL)
 - ligand.sdf            : Ligand SDF for parameterization
 - mdp/*.mdp             : Minimization, NVT, NPT, production parameter files
 - run_gromacs.sh        : Example command sequence (bash)

Suggested toolchain: GROMACS 2021+ with Amber14SB and TIP3P. Ligand via ACPYPE (GAFF2, AM1-BCC).
Note: Keep Mn2+; default ion types in GROMACS cover Mn2+.

Steps:
1) Generate ligand GROMACS parameters with ACPYPE (AmberTools installed)
   acpype -i ligand.sdf -b MOL -a gaff2 -c bcc
   # This creates MOL_GMX.itp, MOL_GMX.gro (or .pdb), etc., in MOL.acpype folder.

2) Build protein topology with pdb2gmx (Amber14SB)
   gmx pdb2gmx -f complex.pdb -o protein.gro -p topol.top -i posre.itp -ff amber14sb -water tip3p
   # Choose Amber14SB; ensure Mn2+ is recognized as MN

3) Merge ligand coordinates into protein.gro (or use complex.pdb as starting coords after conversion):
   # Simplest: convert complex.pdb to GRO
   gmx editconf -f complex.pdb -o complex.gro

4) Define box and solvate
   gmx editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt cubic
   gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top

5) Add ions (0.15 M)
   gmx grompp -f mdp/minim.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
   gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 << EOF
13
EOF

6) Include ligand in topology
   # Edit topol.top: add after forcefield includes
   #   #include "MOL_GMX.itp"
   # In [ system ]/[ molecules ] section, add a line for MOL 1
   # Make sure protein is named appropriately (e.g., Protein_chain_A 1)

7) NVT equilibration
   gmx grompp -f mdp/nvt.mdp -c solv_ions.gro -p topol.top -o nvt.tpr -r solv_ions.gro -maxwarn 2
   gmx mdrun -deffnm nvt

8) NPT equilibration
   gmx grompp -f mdp/npt.mdp -c nvt.gro -p topol.top -o npt.tpr -r nvt.gro -maxwarn 2
   gmx mdrun -deffnm npt

9) Production (10 ns)
   gmx grompp -f mdp/md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 2
   gmx mdrun -deffnm md

Optional: Create index groups for tc-grps matching mdp (Protein_LIG and Water_and_ions):
   gmx make_ndx -f npt.gro -o index.ndx
   # Create group combining Protein and MOL as Protein_LIG, and SOL+NA+CL+MN as Water_and_ions
   # Then set tc-grps in mdp or pass -n index.ndx to grompp

Metal coordination restraints (optional but recommended initially):
 - After building the topology, identify Mn (MN) atom index and coordinating ligand atoms/residues
 - Add a [ distance_restraints ] section in an included .itp with target distances ~2.1–2.3 Å for O–Mn contacts
 - Include that .itp from topol.top

Notes:
 - If ACPYPE outputs different ligand residue name, keep it consistent in .itp and [ molecules ].
 - If pdb2gmx renames your protein, adjust [ molecules ] entries accordingly.
 - Validate the system with gmx check and a short 100 ps test before the 10 ns production.
"""

run_sh = """#!/usr/bin/env bash
set -euo pipefail

# 1) Ligand parameters via ACPYPE (requires AmberTools)
acpype -i ligand.sdf -b MOL -a gaff2 -c bcc || { echo "ACPYPE failed"; exit 1; }

# 2) Protein topology (Amber14SB)
gmx pdb2gmx -f complex.pdb -o protein.gro -p topol.top -i posre.itp -ff amber14sb -water tip3p << EOF
1
EOF

# 3) Convert complex coords to GRO
gmx editconf -f complex.pdb -o complex.gro

# 4) Box and solvate
gmx editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt cubic
gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top

# 5) Ions (0.15 M)
gmx grompp -f mdp/minim.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# 6) Topology edits for ligand (manual step)
echo "Edit topol.top to #include MOL_GMX.itp and add 'MOL 1' under [ molecules ]."

echo "Proceed with NVT/NPT/MD after topology edit."

# 7) NVT
gmx grompp -f mdp/nvt.mdp -c solv_ions.gro -p topol.top -o nvt.tpr -r solv_ions.gro -maxwarn 2
gmx mdrun -deffnm nvt

# 8) NPT
gmx grompp -f mdp/npt.mdp -c nvt.gro -p topol.top -o npt.tpr -r nvt.gro -maxwarn 2
gmx mdrun -deffnm npt

# 9) Production
gmx grompp -f mdp/md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 2
gmx mdrun -deffnm md
"""

# Prepare top-N bundles
OUTROOT.mkdir(exist_ok=True)
sel = pd.read_csv(REDOCK_DETAILED).sort_values('score_ANDV_metal').head(TOP_N)
for _, row in sel.iterrows():
    name = row['Name']
    lig_sdf = DOCK_DIR / row['Path']
    src_complex = MDROOT / name / 'complex_MOL.pdb'
    if not src_complex.exists():
        print(f"Missing complex for {name}, skipping")
        continue
    outdir = OUTROOT / name
    (outdir / 'mdp').mkdir(parents=True, exist_ok=True)
    # Copy files
    shutil.copy2(src_complex, outdir / 'complex.pdb')
    shutil.copy2(lig_sdf, outdir / 'ligand.sdf')
    # Write mdp files
    (outdir / 'mdp' / 'minim.mdp').write_text(minim_mdp)
    (outdir / 'mdp' / 'nvt.mdp').write_text(nvt_mdp)
    (outdir / 'mdp' / 'npt.mdp').write_text(npt_mdp)
    (outdir / 'mdp' / 'md.mdp').write_text(md_mdp)
    # README and script
    (outdir / 'README_GROMACS.txt').write_text(readme_tpl)
    sh = outdir / 'run_gromacs.sh'
    sh.write_text(run_sh)
    os.chmod(sh, 0o755)
    print(f"Prepared GROMACS bundle: {outdir}")

print("All done. See gmx_prep/ for per-complex folders.")
