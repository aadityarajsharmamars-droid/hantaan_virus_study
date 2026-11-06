import os, time, sys
import pandas as pd
from pathlib import Path
from playmolecule.apps.parameterize.v2 import main_parameterize

TOP_N = 10
BASE = Path('.')
REDOCK_DETAILED = BASE / 'potential drugs' / 'redock_metal_modified_scores_detailed.csv'
LIG_DIR = BASE / 'docking' / 'ANDV_6Q99_pocket1_metal_modified'
OUTROOT = BASE / 'md_runs'

assert REDOCK_DETAILED.exists(), f"Missing {REDOCK_DETAILED}"
assert LIG_DIR.exists(), f"Missing {LIG_DIR}"

OUTROOT.mkdir(exist_ok=True)

# Select top N analogues by most negative score
p = pd.read_csv(REDOCK_DETAILED)
p = p.sort_values('score_ANDV_metal', ascending=True).head(TOP_N).copy()

selections = []
for _, row in p.iterrows():
    name = str(row['Name'])  # e.g., outlig57
    sdf_rel = str(row['Path'])  # e.g., outlig57.sdf
    lig_path = LIG_DIR / sdf_rel
    if not lig_path.exists():
        raise FileNotFoundError(f"Ligand SDF not found: {lig_path}")
    selections.append((name, lig_path))

submitted = []
for name, lig_path in selections:
    outdir = OUTROOT / name / 'parameterize'
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Submitting parameterization for {name}: {lig_path}")
    tries = 0
    while True:
        tries += 1
        try:
            job = main_parameterize(outdir=str(outdir), molfile=str(lig_path), forcefield='GAFF2', minimizer='ffeval', charge_type='am1-bcc')
            job.run()
            submitted.append((name, outdir))
            print(f"Submitted {name}")
            break
        except Exception as e:
            if tries >= 5:
                print(f"ERROR: Failed to submit {name} after {tries} attempts: {e}")
                break
            wait = 10 * tries
            print(f"Retry {tries} for {name} after error: {e}\nSleeping {wait}s...", file=sys.stderr)
            time.sleep(wait)

print("Done submitting parameterization jobs. You will receive notifications as each completes.")
for name, outdir in submitted:
    print(f" - {name} -> {outdir}")
