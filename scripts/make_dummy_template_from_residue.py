import sys, os
from statistics import mean

pdb_path = sys.argv[1]
resname = sys.argv[2]
chain = sys.argv[3] if len(sys.argv) > 3 else None
resid = sys.argv[4] if len(sys.argv) > 4 else None
out_sdf = sys.argv[5] if len(sys.argv) > 5 else 'templates/template_center.sdf'

xs, ys, zs = [], [], []
with open(pdb_path, 'r') as f:
    for line in f:
        if not (line.startswith('HETATM') or line.startswith('ATOM')):
            continue
        rn = line[17:20].strip()
        ch = line[21].strip()
        ri = line[22:26].strip()
        if rn != resname:
            continue
        if chain and ch != chain:
            continue
        if resid and ri != resid:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        xs.append(x); ys.append(y); zs.append(z)

if not xs:
    raise SystemExit(f'No atoms found for {resname} {chain or "*"}:{resid or "*"} in {pdb_path}')

cx, cy, cz = mean(xs), mean(ys), mean(zs)
os.makedirs(os.path.dirname(out_sdf), exist_ok=True)
with open(out_sdf, 'w') as w:
    w.write('TEMPLATE_CENTER\n')
    w.write('     RDKit          3D\n\n')
    w.write('  1  0  0  0  0  0  0  0  0  0999 V2000\n')
    w.write(f'{cx:10.4f}{cy:10.4f}{cz:10.4f} C   0  0  0  0  0  0  0  0  0  0  0  0\n')
    w.write('M  END\n$$$$\n')
print(f'Wrote {out_sdf} at ({cx:.3f},{cy:.3f},{cz:.3f})')
