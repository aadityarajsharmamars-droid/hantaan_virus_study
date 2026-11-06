import os
import sys
import requests
from pathlib import Path
from typing import List, Tuple

UA = {'User-Agent': 'PlayMoleculeAI/phylogeny (contact: user)'}

# Targets: (label, UniProt query)
TARGETS = [
    ("Influenza_A_PB1", 'taxonomy_id:11320 AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Influenza_B_PB1", 'taxonomy_id:11520 AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Influenza_C_PB1", 'taxonomy_id:11552 AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Andes_orthohantavirus_L", 'organism_name:"Andes orthohantavirus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Hantaan_virus_L", 'organism_name:"Hantaan virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Puumala_virus_L", 'organism_name:"Puumala virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Rabies_lyssavirus_L", 'organism_name:"Rabies virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)')
]

BASE = 'https://rest.uniprot.org/uniprotkb/search'


def fetch_one_fasta(query: str, reviewed_first: bool = True) -> Tuple[str, str]:
    """Return (header, sequence) for first hit in FASTA format."""
    # Try reviewed first if requested
    queries = []
    if reviewed_first:
        queries.append(f"({query}) AND reviewed:true")
    queries.append(query)  # fallback unfiltered

    for q in queries:
        params = {
            'query': q,
            'format': 'fasta',
            'size': 1
        }
        r = requests.get(BASE, params=params, headers=UA, timeout=60)
        if r.status_code == 200 and r.text.strip():
            fasta = r.text.strip().splitlines()
            if fasta and fasta[0].startswith('>'):
                header = fasta[0]
                seq = ''.join(line.strip() for line in fasta[1:] if line and not line.startswith('>'))
                if len(seq) >= 200:  # ensure polymerase-length; PB1 ~ 757 aa, L > 2000 aa; allow shorter for PA? but here PB1/L
                    return header, seq
                else:
                    # might have retrieved short fragments; continue
                    continue
    raise RuntimeError(f"No sequence found for query: {query}")


def main(outdir: str = 'phylogeny_rdRp'):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    seq_records = []  # list of (label, uniprot_header, sequence)

    for label, q in TARGETS:
        try:
            header, seq = fetch_one_fasta(q, reviewed_first=True)
            seq_records.append((label, header, seq))
            print(f"Fetched {label}: {header}")
        except Exception as e:
            print(f"WARNING: failed to fetch {label} with query '{q}': {e}")

    if len(seq_records) < 4:
        print("ERROR: insufficient sequences fetched for tree construction.")
        sys.exit(1)

    # Write FASTA with clean labels to avoid long names
    fasta_path = out / 'rdRp_sequences.fasta'
    with fasta_path.open('w') as fh:
        for label, header, seq in seq_records:
            fh.write(f">{label}\n")
            # wrap sequence
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
    print(f"Wrote {fasta_path}")

    # Align with MAFFT
    aligned_path = out / 'rdRp_aligned.fasta'
    cmd = f"mafft --auto --thread -1 {fasta_path} > {aligned_path} 2> {out/'mafft.log'}"
    ret = os.system(cmd)
    if ret != 0 or not aligned_path.exists():
        print("ERROR: MAFFT alignment failed.")
        sys.exit(2)
    print(f"Aligned -> {aligned_path}")

    # Build NJ tree using Biopython
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

    aln = AlignIO.read(str(aligned_path), 'fasta')
    calc = DistanceCalculator('blosum62')
    dm = calc.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Root with outgroup if present
    outgroup = None
    for clade in tree.find_clades():
        if clade.name == 'Rabies_lyssavirus_L':
            outgroup = clade
            break
    if outgroup is not None:
        tree.root_with_outgroup(outgroup)

    newick_path = out / 'rdRp_tree.nwk'
    Phylo.write(tree, str(newick_path), 'newick')
    print(f"Tree -> {newick_path}")

    # Render a simple figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8, 6), dpi=150)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    fig.tight_layout()
    fig_path = out / 'rdRp_tree.png'
    fig.savefig(fig_path)
    print(f"Figure -> {fig_path}")


if __name__ == '__main__':
    outdir = sys.argv[1] if len(sys.argv) > 1 else 'phylogeny_rdRp'
    main(outdir)
