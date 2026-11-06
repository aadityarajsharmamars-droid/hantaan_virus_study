import os
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

BASE = Path('.')
RECEPTOR = BASE / 'receptors' / '6Q99' / 'output_metal.pdb'
REDOCK_DETAILED = BASE / 'potential drugs' / 'redock_metal_modified_scores_detailed.csv'
DOCK_DIR = BASE / 'docking' / 'ANDV_6Q99_pocket1_metal_modified'
OUTROOT = BASE / 'md_runs'
TOP_N = 10

assert RECEPTOR.exists(), f'Missing receptor: {RECEPTOR}'
assert REDOCK_DETAILED.exists()

# Load top N selections
p = pd.read_csv(REDOCK_DETAILED).sort_values('score_ANDV_metal').head(TOP_N)
selections = [(row['Name'], DOCK_DIR / row['Path']) for _, row in p.iterrows()]

for name, sdf_path in selections:
    outdir = OUTROOT / name
    outdir.mkdir(parents=True, exist_ok=True)
    complex_pdb = outdir / 'complex_MOL.pdb'
    lig_pdb = outdir / 'ligand_MOL.pdb'

    # Load ligand first conformer
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    if len(suppl) == 0 or suppl[0] is None:
        raise ValueError(f'Cannot read SDF: {sdf_path}')
    mol = Chem.Mol(suppl[0])
    # Assign PDB residue info with resname MOL
    for i, atom in enumerate(mol.GetAtoms()):
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName('MOL')
        info.SetResidueNumber(500)
        info.SetChainId('L')
        info.SetName((atom.GetSymbol() + '  ')[:4])
        atom.SetPDBResidueInfo(info)
    pdb_block = Chem.MolToPDBBlock(mol)
    lig_pdb.write_text(pdb_block)

    # Merge receptor and ligand PDBs by concatenation
    rec_txt = RECEPTOR.read_text()
    # Ensure receptor ends with ENDMDL/TER if needed
    if not rec_txt.strip().endswith('ENDMDL'):
        rec_txt = rec_txt.rstrip() + '\nENDMDL\n'
    with open(complex_pdb, 'w') as f:
        f.write(rec_txt)
        # Remove headers from ligand block except ATOM/HETATM lines
        lig_lines = [ln for ln in pdb_block.splitlines() if ln.startswith(('ATOM', 'HETATM'))]
        # Convert ATOM to HETATM for ligand
        lig_lines = [ ('HETATM' + ln[6:]) if ln.startswith('ATOM') else ln for ln in lig_lines ]
        for ln in lig_lines:
            f.write(ln + '\n')
        f.write('END\n')
    print(f'Wrote complex: {complex_pdb}')
