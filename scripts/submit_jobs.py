import time
from playmolecule.apps.proteinprepare import proteinprepare
from playmolecule.apps.deepsite import deepsite

jobs = [
    ('proteinprepare', {'outdir':'receptors/6FS6','pdbid':'6FS6','pH':7.2}),
    ('proteinprepare', {'outdir':'receptors/6Q99','pdbid':'6Q99','pH':7.2}),
    ('deepsite', {'outdir':'pockets/6FS6','pdbid':'6FS6'}),
    ('deepsite', {'outdir':'pockets/6Q99','pdbid':'6Q99'}),
]

for name, kwargs in jobs:
    for attempt in range(5):
        try:
            if name=='proteinprepare':
                job = proteinprepare(**kwargs)
            else:
                job = deepsite(**kwargs)
            job.run()
            print(f"Submitted {name} with {kwargs}")
            break
        except Exception as e:
            wait = 5*(attempt+1)
            print(f"Error submitting {name} {kwargs}: {e}. Retrying in {wait}s")
            time.sleep(wait)
    time.sleep(5)
print('Done submitting.')
