from playmolecule.apps.proteinprepare import proteinprepare
from playmolecule.apps.deepsite import deepsite
import os

os.makedirs('receptors/6FS6', exist_ok=True)
os.makedirs('receptors/6Q99', exist_ok=True)
os.makedirs('pockets/6FS6', exist_ok=True)
os.makedirs('pockets/6Q99', exist_ok=True)

# Submit protein preparation
job1 = proteinprepare(outdir='receptors/6FS6', pdbid='6FS6', pH=7.2)
job1.run()
job2 = proteinprepare(outdir='receptors/6Q99', pdbid='6Q99', pH=7.2)
job2.run()

# Submit pocket detection using PDB IDs (independent of ProteinPrepare outputs)
job3 = deepsite(outdir='pockets/6FS6', pdbid='6FS6')
job3.run()
job4 = deepsite(outdir='pockets/6Q99', pdbid='6Q99')
job4.run()

print('Submitted ProteinPrepare for 6FS6 and 6Q99, and DeepSite for both.')
