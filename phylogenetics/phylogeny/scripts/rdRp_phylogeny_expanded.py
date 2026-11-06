import os
import re
import sys
import requests
from pathlib import Path
from typing import List, Tuple

UA = {'User-Agent': 'PlayMoleculeAI/phylogeny (contact: user)'}
BASE = 'https://rest.uniprot.org/uniprotkb/search'

# Expanded representative list across Articulavirales and Bunyavirales, plus Mononegavirales outgroups
TARGETS = [
    ("Influenza_A_PB1", 'taxonomy_id:11320 AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Influenza_B_PB1", 'taxonomy_id:11520 AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Influenza_C_PB1", 'taxonomy_id:11552 AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Influenza_D_PB1", 'organism_name:"Influenza D virus" AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Isavirus_PB1", 'organism_name:"Infectious salmon anemia virus" AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Thogoto_virus_PB1", 'organism_name:"Thogoto virus" AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),
    ("Quaranfil_virus_PB1", 'organism_name:"Quaranfil virus" AND (gene_exact:PB1 OR protein_name:"RNA-directed RNA polymerase catalytic subunit")'),

    ("Andes_orthohantavirus_L", 'organism_name:"Andes orthohantavirus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Hantaan_virus_L", 'organism_name:"Hantaan virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Puumala_virus_L", 'organism_name:"Puumala virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Bunyamwera_virus_L", 'organism_name:"Bunyamwera virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("La_Crosse_virus_L", 'organism_name:"La Crosse virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Rift_Valley_fever_virus_L", 'organism_name:"Rift Valley fever virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Crimean_Congo_hemorrhagic_fever_virus_L", 'organism_name:"Crimean-Congo hemorrhagic fever virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Tomato_spotted_wilt_virus_L", 'organism_name:"Tomato spotted wilt virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),

    ("Rabies_lyssavirus_L", 'organism_name:"Rabies virus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
    ("Measles_morbillivirus_L", 'organism_name:"Measles morbillivirus" AND (protein_name:polymerase OR protein_name:L OR gene_exact:L)'),
]


def fetch_one_fasta(query: str, reviewed_first: bool = True) -> Tuple[str, str]:
    queries = []
    if reviewed_first:
        queries.append(f"({query}) AND reviewed:true")
    queries.append(query)

    for q in queries:
        params = {'query': q, 'format': 'fasta', 'size': 1}
        r = requests.get(BASE, params=params, headers=UA, timeout=60)
        if r.status_code == 200 and r.text.strip():
            fasta = r.text.strip().splitlines()
            if fasta and fasta[0].startswith('>'):
                header = fasta[0]
                seq = ''.join(line.strip() for line in fasta[1:] if line and not line.startswith('>'))
                if len(seq) >= 200:
                    return header, seq
    raise RuntimeError(f"No sequence found for query: {query}")


def parse_year_from_header(header: str):
    # look for 4-digit years 1900-2025
    yrs = re.findall(r'(19[0-9]{2}|20[0-2][0-9])', header)
    if yrs:
        # choose the oldest plausible sampling year in header
        return int(sorted(yrs)[0])
    return None


def main(outdir: str = 'phylogeny_expanded'):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    meta = []

    seq_records = []
    for label, q in TARGETS:
        try:
            header, seq = fetch_one_fasta(q, reviewed_first=True)
            year = parse_year_from_header(header)
            seq_records.append((label, header, seq))
            meta.append((label, header, year))
            print(f"Fetched {label}: year={year} | {header}")
        except Exception as e:
            print(f"WARNING: failed to fetch {label} with query '{q}': {e}")

    if len(seq_records) < 8:
        print("ERROR: insufficient sequences fetched.")
        sys.exit(1)

    # Write FASTA
    fasta_path = out / 'rdRp_sequences_expanded.fasta'
    with fasta_path.open('w') as fh:
        for label, header, seq in seq_records:
            fh.write(f">{label}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
    print(f"Wrote {fasta_path}")

    # Save metadata TSV
    with (out / 'metadata.tsv').open('w') as fh:
        fh.write('label\theader\tyear\n')
        for label, header, year in meta:
            fh.write(f"{label}\t{header}\t{year if year is not None else ''}\n")

    # Align with MAFFT
    aligned_path = out / 'rdRp_aligned_expanded.fasta'
    cmd = f"mafft --auto --thread -1 {fasta_path} > {aligned_path} 2> {out/'mafft.log'}"
    ret = os.system(cmd)
    if ret != 0 or not aligned_path.exists():
        print("ERROR: MAFFT alignment failed.")
        sys.exit(2)
    print(f"Aligned -> {aligned_path}")

    # Build NJ tree
    from Bio import AlignIO, Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

    aln = AlignIO.read(str(aligned_path), 'fasta')
    calc = DistanceCalculator('blosum62')
    dm = calc.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Root with outgroups if present (Mononegavirales)
    outgroups = [c for c in tree.find_clades() if getattr(c, 'name', '') in ('Rabies_lyssavirus_L','Measles_morbillivirus_L')]
    if outgroups:
        tree.root_with_outgroup(outgroups[0])

    newick_path = out / 'rdRp_tree_expanded.nwk'
    Phylo.write(tree, str(newick_path), 'newick')
    print(f"Tree -> {newick_path}")

    # Draw figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(9, 7), dpi=150)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    fig.tight_layout()
    fig.savefig(out / 'rdRp_tree_expanded.png')

    # Attempt crude root-to-tip regression if enough dated tips
    dated = [(label, year) for label, _, year in meta if year is not None]
    if len(dated) >= 4:
        # compute path lengths from root to each dated tip
        root = tree.root
        # map tip names to clades
        tips = {clade.name: clade for clade in tree.get_terminals()}
        x = []  # years
        y = []  # root-to-tip distances
        used = []
        for label, year in dated:
            clade = tips.get(label)
            if clade is None:
                continue
            # compute distance
            dist = tree.distance(root, clade)
            if dist is not None:
                x.append(year)
                y.append(dist)
                used.append(label)
        if len(x) >= 4:
            import numpy as np
            X = np.vstack([np.array(x), np.ones(len(x))]).T
            beta, intercept = np.linalg.lstsq(X, np.array(y), rcond=None)[0]
            rate = beta  # substitutions per unit year (arbitrary units from NJ distances)
            t_root = -intercept / beta if beta != 0 else None
            with (out / 'clock_regression.txt').open('w') as fh:
                fh.write(f"Used tips: {used}\n")
                fh.write(f"Slope (rate, arb. units/year): {rate}\n")
                fh.write(f"Intercept: {intercept}\n")
                fh.write(f"Estimated root time (year): {t_root}\n")
            # quick plot
            import matplotlib.pyplot as plt
            plt.figure(figsize=(5,4), dpi=150)
            plt.scatter(x,y)
            if beta!=0:
                xs = np.linspace(min(x)-10, max(x)+10, 100)
                ys = beta*xs + intercept
                plt.plot(xs, ys, 'r--')
            plt.xlabel('Sampling year')
            plt.ylabel('Root-to-tip distance (NJ)')
            plt.tight_layout()
            plt.savefig(out / 'clock_regression.png')
            print('Clock regression written.')
        else:
            print('Not enough dated tips for regression after mapping.')
    else:
        print('Insufficient tip dates for time calibration. See metadata.tsv and consider external BEAST/TreeTime analysis.')

if __name__ == '__main__':
    outdir = sys.argv[1] if len(sys.argv) > 1 else 'phylogeny_expanded'
    main(outdir)
