import time
from playmolecule.apps import acedock

jobs = []

# Dock to Influenza PA endonuclease (6FS6) at BXA center
jobs.append((
    'docking/PA_6FS6',
    {
        'outdir': 'docking/PA_6FS6',
        'protein': 'receptors/6FS6/output_noE4Z.pdb',
        'template': 'templates/6FS6_bxa_center.sdf',
        'ligands': 'ligands/ligands.sdf',
        'mode': 'free',
        'radius': 8,
        'n_poses': 10,
        'rescore_pharma': False,
        'dummy_template': True
    }
))

# Dock to Andes L endonuclease pocket 1 center
jobs.append((
    'docking/ANDV_6Q99_pocket1',
    {
        'outdir': 'docking/ANDV_6Q99_pocket1',
        'protein': 'receptors/6Q99/output.pdb',
        'template': 'dock_inputs/pocket1_center_template.sdf',
        'ligands': 'ligands/ligands.sdf',
        'mode': 'free',
        'radius': 8,
        'n_poses': 10,
        'rescore_pharma': False,
        'dummy_template': True
    }
))

for outdir, kwargs in jobs:
    for attempt in range(4):
        try:
            job = acedock.v1.acedock(**kwargs)
            job.run()
            print(f"Submitted AceDock job: {outdir}")
            break
        except Exception as e:
            wait = 10*(attempt+1)
            print(f"Error submitting {outdir}: {e}. Retrying in {wait}s")
            time.sleep(wait)
    time.sleep(5)
print('All docking jobs submitted.')
