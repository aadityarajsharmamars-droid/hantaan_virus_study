import sys
inp = sys.argv[1]
outp = sys.argv[2]
resname = sys.argv[3]
chain = sys.argv[4] if len(sys.argv) > 4 else None
resid = sys.argv[5] if len(sys.argv) > 5 else None

with open(inp,'r') as f, open(outp,'w') as w:
    for line in f:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            rn = line[17:20].strip()
            ch = line[21].strip()
            ri = line[22:26].strip()
            if rn == resname and (chain is None or ch==chain) and (resid is None or ri==resid):
                continue
        w.write(line)
print(f'Wrote {outp}')
