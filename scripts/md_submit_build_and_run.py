import os, sys, time, zipfile
from pathlib import Path
import pandas as pd
from playmolecule.apps.systembuilder import systembuilder
from playmolecule.apps.acemd import simplerun

BASE = Path('.')
OUTROOT = BASE / 'md_runs'
TOP_N = 10

# Helper: robust submission

def submit_with_retry(make_job_func, label, max_tries=5):
    tries = 0
    while True:
        tries += 1
        try:
            job = make_job_func()
            job.run()
            print(f"Submitted {label}")
            return True
        except Exception as e:
            if tries >= max_tries:
                print(f"ERROR: Failed to submit {label} after {tries} attempts: {e}")
                return False
            wait = 10 * tries
            print(f"Retry {tries} for {label} after error: {e}\nSleeping {wait}s...", file=sys.stderr)
            time.sleep(wait)

submitted = []
for ligdir in sorted([p for p in OUTROOT.glob('outlig*') if p.is_dir()])[:TOP_N]:
    name = ligdir.name
    print(f"Processing {name}")
    complex_pdb = ligdir / 'complex_MOL.pdb'
    param_dir = ligdir / 'parameterize' / 'parameters'
    if not complex_pdb.exists():
        print(f" - Missing complex PDB: {complex_pdb}, skipping")
        continue
    if not param_dir.exists():
        print(f" - Parameters not ready yet: {param_dir}, skipping")
        continue
    # Find parameter files
    frcmods = list(param_dir.glob('*.frcmod'))
    cifs = [p for p in param_dir.glob('*-orig.cif')] + [p for p in param_dir.glob('*.cif')]
    if not frcmods or not cifs:
        print(f" - Incomplete parameters in {param_dir}, skipping")
        continue
    frcmod = frcmods[0]
    cif = cifs[0]
    # Prepare zip with MOL filenames
    ligzip = ligdir / 'MOL_params.zip'
    with zipfile.ZipFile(ligzip, 'w', zipfile.ZIP_DEFLATED) as z:
        z.write(frcmod, arcname='MOL.frcmod')
        z.write(cif, arcname='MOL.cif')
    print(f" - Packed params: {ligzip}")

    # Submit SystemBuilder
    build_out = ligdir / 'build'
    build_out.mkdir(exist_ok=True)
    def make_build_job():
        return systembuilder(outdir=str(build_out), pdb=str(complex_pdb), pH=7.4, ligands=[str(ligzip)], saltconc=0.15)
    ok = submit_with_retry(make_build_job, f"SystemBuilder {name}")
    if not ok:
        continue

    # Submit two replicas of SimpleRun (they will queue and start after build completes)
    for rep in (1, 2):
        sim_out = ligdir / f'sim_r{rep}'
        sim_out.mkdir(exist_ok=True)
        def make_sim_job(rep=rep):
            return simplerun(outdir=str(sim_out), builddir=str(build_out), equilrun=1.0, prodrun=10.0, temperature=300.0, ligresnames=['MOL'], defaultrestraints=True, restraintdecay=0.5)
        ok2 = submit_with_retry(make_sim_job, f"SimpleRun {name} rep{rep}")
        if ok2:
            submitted.append((name, rep))

print("Submission pipeline finished.")
for name, rep in submitted:
    print(f" - {name} sim_r{rep}")
