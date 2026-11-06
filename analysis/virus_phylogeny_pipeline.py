import os
import time
import csv
import re
from pathlib import Path
from typing import List, Dict, Tuple
from urllib.parse import quote

from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import Consensus

# Configure Entrez (NCBI E-utilities)
Entrez.email = "playmoleculeai@acellera.com"
Entrez.tool = "PlayMoleculeAI"

PROJECT = Path('.')
RAW_DIR = PROJECT / 'proteins_raw'
OUT_DIR = PROJECT / 'domain_phylogeny'
RAW_DIR.mkdir(exist_ok=True, parents=True)
OUT_DIR.mkdir(exist_ok=True, parents=True)

VIRUSES = [
    'Andes virus', 'Hantaan virus', 'Seoul virus', 'Sin Nombre virus', 'Puumala virus', 'Dobrava-Belgrade virus', 'Tula virus', 'Prospect Hill virus',
    'La Crosse virus', 'California encephalitis virus', 'Bunyamwera virus', 'Oropouche virus', 'Schmallenberg virus', 'Akabane virus',
    'Rift Valley fever virus', 'SFTS virus', 'Heartland virus', 'Uukuniemi virus', 'Toscana virus', 'Crimean-Congo hemorrhagic fever virus',
    'Nairobi sheep disease virus', 'Dugbe virus', 'Hazara virus',
    'Lassa virus', 'LCMV', 'Junin virus', 'Machupo virus', 'Guanarito virus',
    'Influenza A virus', 'Influenza B virus', 'Influenza C virus', 'Influenza D virus', 'Thogoto virus', 'Dhori virus', 'Quaranfil virus', 'Infectious salmon anemia virus',
    'Vesicular stomatitis virus', 'Rabies virus', 'Measles virus', 'Ebola virus'
]

# Map each virus to the desired proteins for download
# Orthomyxoviridae: PA (endonuclease) + PB1 (polymerase)
ORTHO = {'Influenza A virus', 'Influenza B virus', 'Influenza C virus', 'Influenza D virus', 'Thogoto virus', 'Dhori virus', 'Quaranfil virus', 'Infectious salmon anemia virus'}

# Rhabdoviridae/Paramyxoviridae/Filoviridae: L
L_ONLY = {'Vesicular stomatitis virus', 'Rabies virus', 'Measles virus', 'Ebola virus'}

# Arenaviridae + Bunyavirales: L
# Explicit list for clarity
L_LIST = set(VIRUSES) - ORTHO - L_ONLY

# Build per-virus target protein names
TARGETS: Dict[str, List[str]] = {}
for v in VIRUSES:
    if v in ORTHO:
        TARGETS[v] = ['PA', 'PB1']
    else:
        TARGETS[v] = ['L']


def esearch(db: str, term: str, retmax: int = 10):
    for attempt in range(3):
        try:
            with Entrez.esearch(db=db, term=term, retmax=retmax) as h:
                out = Entrez.read(h)
                time.sleep(0.3)
                return out
        except Exception:
            time.sleep(1 + attempt)
    raise RuntimeError(f"esearch failed for {db}: {term}")


def efetch_protein(accession: str, fmt: str = 'fasta') -> str:
    for attempt in range(3):
        try:
            with Entrez.efetch(db='protein', id=accession, rettype=fmt, retmode='text') as h:
                return h.read()
        except Exception as e:
            time.sleep(1 + attempt)
    raise RuntimeError(f"efetch protein failed: {accession}")


def fetch_one_protein(virus: str, target: str) -> Tuple[str, str, int, str]:
    # Build query per target
    if target == 'PA':
        key = '(PA[Gene] OR polymerase acidic[Title] OR endonuclease[Title])'
        minlen = 180
    elif target == 'PB1':
        key = '(PB1[Gene] OR RNA-directed RNA polymerase catalytic subunit[Title] OR polymerase[Title])'
        minlen = 600
    else:  # L
        key = '(L[Gene] OR RNA-dependent RNA polymerase[Title] OR polymerase[Title])'
        minlen = 1200
    q = f'({virus}[Organism]) AND {key}'
    res = esearch('protein', q, retmax=50)
    ids = res.get('IdList', [])
    acc = None
    seq = None
    source = 'NCBI Protein'
    from io import StringIO
    for uid in ids:
        try:
            fa = efetch_protein(uid, 'fasta')
        except Exception:
            continue
        if not fa or not fa.startswith('>'):
            continue
        handle = StringIO(fa)
        records = list(SeqIO.parse(handle, 'fasta'))
        if not records:
            continue
        r = records[0]
        if len(r.seq) < minlen:
            continue
        header = fa.splitlines()[0]
        if '|' in header:
            parts = header.split('|')
            acc = parts[1]
        else:
            acc = r.id
        seq = fa
        break
    if acc is None or seq is None:
        raise RuntimeError(f"No protein found for {virus} {target}")
    length = len(records[0].seq)
    return acc, seq, length, source


def fetch_taxonomy(virus: str) -> Tuple[str, Dict[str, str]]:
    # Resolve taxid
    sr = esearch('taxonomy', virus, retmax=1)
    if not sr.get('IdList'):
        return '', {}
    taxid = sr['IdList'][0]
    for attempt in range(3):
        try:
            with Entrez.efetch(db='taxonomy', id=taxid) as h:
                rec = Entrez.read(h)[0]
                lineage = rec.get('LineageEx', [])
                ranks = {x['Rank']: x['ScientificName'] for x in lineage if x.get('Rank')}
                # add self at species
                ranks.setdefault('species', rec.get('ScientificName',''))
                return taxid, ranks
        except Exception:
            time.sleep(1 + attempt)
    return taxid, {}


def save_metadata(rows: List[Dict[str, str]]):
    out_csv = OUT_DIR / 'virus_metadata.csv'
    cols = ['Virus_name','Family','Order','Class','Phylum','NCBI_taxid','Protein_name','Accession_ID','Protein_length','Downloaded_from']
    with out_csv.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, '') for c in cols})
    print('Metadata ->', out_csv)


def extract_endonuclease_domains(raw_fasta_paths: List[Path]) -> Path:
    # Heuristic: for PA or L, take N-terminal 200 aa as endonuclease
    out_fa = OUT_DIR / 'endonuclease_domains.fasta'
    with out_fa.open('w') as out:
        for fp in raw_fasta_paths:
            for rec in SeqIO.parse(str(fp), 'fasta'):
                name = fp.stem
                seq = str(rec.seq)
                # determine if PA or L
                prot = 'PA' if 'PA' in name.split('__')[-1] else 'L'
                dom = seq[:200]
                if len(dom) < 100:
                    continue
                out.write(f">{name}\n")
                for i in range(0, len(dom), 60):
                    out.write(dom[i:i+60]+'\n')
    print('Endonuclease domains ->', out_fa)
    return out_fa


def align_and_trim(in_fa: Path, prefix: str) -> Tuple[Path, Path]:
    aln_fa = OUT_DIR / f'{prefix}_aligned.fasta'
    trimmed_fa = OUT_DIR / f'{prefix}_trimmed.fasta'
    os.system(f"mafft --maxiterate 1000 --localpair {in_fa} > {aln_fa} 2> {OUT_DIR/'mafft.log'}")
    # simple trim: remove columns with >50% gaps
    aln = AlignIO.read(str(aln_fa), 'fasta')
    keep_cols = []
    for i in range(aln.get_alignment_length()):
        col = aln[:, i]
        gaps = col.count('-') + col.count('.')
        if gaps/len(aln) <= 0.5:
            keep_cols.append(i)
    trimmed = None
    for j, i in enumerate(keep_cols):
        trimmed = aln[:, i:i+1] if j==0 else trimmed + aln[:, i:i+1]
    AlignIO.write(trimmed, str(trimmed_fa), 'fasta')
    print('Trimmed ->', trimmed_fa)
    return aln_fa, trimmed_fa


def build_tree(fa: Path, out_prefix: str, bootstrap: int = 0):
    aln = AlignIO.read(str(fa), 'fasta')
    calc = DistanceCalculator('blosum62')
    constructor = DistanceTreeConstructor(calc, method='nj')
    tree = constructor.build_tree(aln)
    newick = OUT_DIR / f'{out_prefix}.nwk'
    Phylo.write(tree, str(newick), 'newick')
    if bootstrap and len(aln) >= 4:
        cons = Consensus.bootstrap_consensus(aln, bootstrap, constructor, Consensus.majority_consensus)
        Phylo.write(cons, str(OUT_DIR / f'{out_prefix}.bootstrap.nwk'), 'newick')
    # figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10,8), dpi=150)
    ax = fig.add_subplot(1,1,1)
    Phylo.draw(tree, do_show=False, axes=ax)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f'{out_prefix}.png')
    print('Tree ->', newick)


def taxonomy_tree(metadata_csv: Path):
    # Build a simple Newick grouping by Order then Family
    rows = list(csv.DictReader(open(metadata_csv)))
    by_order: Dict[str, Dict[str, List[str]]] = {}
    for r in rows:
        order = r.get('Order','Unknown').replace(' ','_') or 'Unknown'
        fam = r.get('Family','Unknown').replace(' ','_') or 'Unknown'
        sp = r['Virus_name'].replace(' ','_')+"_"+r['Protein_name']
        by_order.setdefault(order, {}).setdefault(fam, []).append(sp)
    def fam_to_newick(d: Dict[str, List[str]]):
        return '('+','.join('('+','.join(v)+')'+f'{k}' for k,v in d.items())+')'
    newick = '(' + ','.join(fam_to_newick(f)+'%s'%ordr for ordr,f in by_order.items()) + ');'
    out = OUT_DIR / 'taxonomy_tree.nwk'
    out.write_text(newick)
    print('Taxonomy tree (grouping) ->', out)


def main():
    metadata_rows = []
    fasta_paths = []

    for virus in VIRUSES:
        taxid, ranks = fetch_taxonomy(virus)
        family = ranks.get('family','')
        order = ranks.get('order','')
        vclass = ranks.get('class','')
        phylum = ranks.get('phylum','')
        for prot in TARGETS[virus]:
            try:
                # Skip if already downloaded for resume capability
                existing = list(RAW_DIR.glob(f"{virus.replace(' ','_')}__{prot}__*.fasta"))
                if existing:
                    fname = existing[0]
                    # estimate length
                    rec = next(SeqIO.parse(str(fname), 'fasta'))
                    length = len(rec.seq)
                    acc = fname.stem.split('__')[-1]
                    source = 'NCBI Protein'
                else:
                    acc, fa, length, source = fetch_one_protein(virus, prot)
                    fname = RAW_DIR / f"{virus.replace(' ','_')}__{prot}__{acc}.fasta"
                    fname.write_text(fa)
                fasta_paths.append(fname)
                metadata_rows.append({
                    'Virus_name': virus,
                    'Family': family,
                    'Order': order,
                    'Class': vclass,
                    'Phylum': phylum,
                    'NCBI_taxid': taxid,
                    'Protein_name': prot,
                    'Accession_ID': acc,
                    'Protein_length': str(length),
                    'Downloaded_from': source
                })
                print(f"OK {virus} {prot} -> {acc} ({length} aa)")
            except Exception as e:
                print(f"FAIL {virus} {prot}: {e}")
    # Save metadata
    save_metadata(metadata_rows)

    # Extract endonuclease domains
    endo_fa = extract_endonuclease_domains(fasta_paths)

    # Align & trim
    aln_fa, trimmed_fa = align_and_trim(endo_fa, prefix='endo')

    # Build trees
    build_tree(trimmed_fa, out_prefix='endo_tree', bootstrap=100)

    # Polymerase tree using full sequences (all PB1 and L)
    poly_fa = OUT_DIR / 'polymerase_full.fasta'
    with poly_fa.open('w') as out:
        for fp in fasta_paths:
            for rec in SeqIO.parse(str(fp), 'fasta'):
                seq = str(rec.seq)
                if len(seq) < 200:
                    continue
                out.write(f">{fp.stem}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60]+'\n')
    _, poly_trim = align_and_trim(poly_fa, prefix='poly')
    build_tree(poly_trim, out_prefix='poly_tree', bootstrap=0)

    # Taxonomy grouping tree
    taxonomy_tree(OUT_DIR / 'virus_metadata.csv')

if __name__ == '__main__':
    main()
