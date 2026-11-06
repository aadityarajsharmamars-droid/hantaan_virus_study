import os
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import Consensus
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base = Path('phylogeny_expanded_v2')
ali_fa = base / 'rdRp_aligned_expanded.fasta'
outdir = Path('phylogeny_bootstrap')
outdir.mkdir(exist_ok=True, parents=True)

# Load alignment
aln = AlignIO.read(str(ali_fa), 'fasta')

# Optionally, filter columns with >70% gaps to reduce noise
keep_cols = []
for i in range(aln.get_alignment_length()):
    col = aln[:, i]
    gaps = col.count('-') + col.count('.')
    if gaps / len(aln) <= 0.7:
        keep_cols.append(i)
if len(keep_cols) < aln.get_alignment_length():
    trimmed = None
    for j, i in enumerate(keep_cols):
        if trimmed is None:
            trimmed = aln[:, i:i+1]
        else:
            trimmed += aln[:, i:i+1]
    aln = trimmed

# Build base NJ tree
calc = DistanceCalculator('blosum62')
constructor = DistanceTreeConstructor(calc, method='nj')
base_tree = constructor.build_tree(aln)

# Bootstrap (500 replicates)
boot_reps = 200
cons_tree = Consensus.bootstrap_consensus(aln, boot_reps, constructor, Consensus.majority_consensus)

# Root with mononegavirales outgroup if present
outgroups = [c for c in cons_tree.find_clades() if getattr(c, 'name','') in ('Rabies_lyssavirus_L','Measles_morbillivirus_L')]
if outgroups:
    cons_tree.root_with_outgroup(outgroups[0])

# Write trees
Phylo.write(base_tree, str(outdir/'tree_nj_no_bootstrap.nwk'), 'newick')
Phylo.write(cons_tree, str(outdir/'tree_nj_bootstrap_consensus.nwk'), 'newick')

# Plot consensus with support values on internal nodes
fig = plt.figure(figsize=(10,8), dpi=160)
ax = fig.add_subplot(1,1,1)
Phylo.draw(cons_tree, do_show=False, axes=ax)
# Add bootstrap values if available
for clade in cons_tree.find_clades():
    if clade.confidence is not None and not clade.is_terminal():
        x = cons_tree.distance(cons_tree.root, clade)
        y = None  # y is not directly accessible; annotate via handle text positions is complex.
        # Instead, set the label attribute so draw function shows it
        clade.name = f"{int(round(clade.confidence))}"
plt.tight_layout()
plt.savefig(outdir/'tree_nj_bootstrap_consensus.png')

print('Bootstrap consensus tree written to', outdir)
