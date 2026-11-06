import os
import re
import csv
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

UA = {'User-Agent': 'PlayMoleculeAI/Negarnaviricota-pipeline (contact: user)'}
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

# Step 1: Virus list (40)
VIRUSES = [
    'Andes virus','Hantaan virus','Seoul virus','Sin Nombre virus','Puumala virus',
    'Dobrava-Belgrade virus','Tula virus','Prospect Hill virus',
    'La Crosse virus','California encephalitis virus','Bunyamwera virus','Oropouche virus',
    'Schmallenberg virus','Akabane virus','Rift Valley fever virus','SFTS virus','Heartland virus',
    'Uukuniemi virus','Toscana virus','Crimean-Congo hemorrhagic fever virus','Nairobi sheep disease virus',
    'Dugbe virus','Hazara virus','Lassa virus','LCMV','Junin virus','Machupo virus','Guanarito virus',
    'Influenza A virus','Influenza B virus','Influenza C virus','Influenza D virus','Thogoto virus','Dhori virus',
    'Quaranfil virus','Infectious salmon anemia virus','Vesicular stomatitis virus','Rabies virus','Measles virus','Ebola virus']

# Aliases to improve taxonomy/protein queries
ALIASES = {
    'SFTS virus': ['severe fever with thrombocytopenia syndrome virus','Huaiyangshan banyangvirus'],
    'LCMV': ['Lymphocytic choriomeningitis mammarenavirus','Lymphocytic choriomeningitis virus'],
    'Infectious salmon anemia virus': ['Infectious salmon anaemia virus','Isavirus'],
    'Rabies virus': ['Rabies lyssavirus'],
    'Measles virus': ['Measles morbillivirus'],
    'Ebola virus': ['Zaire ebolavirus','Ebolavirus'],
}

# Which proteins to retrieve per virus
ORTHOMYXO = set(['Influenza A virus','Influenza B virus','Influenza C virus','Influenza D virus','Thogoto virus','Dhori virus','Quaranfil virus','Infectious salmon anemia virus'])

PROTEIN_TARGETS: Dict[str, List[str]] = {}
for v in VIRUSES:
    if v in ORTHOMYXO:
        PROTEIN_TARGETS[v] = ['PA','PB1']
    else:
        PROTEIN_TARGETS[v] = ['L']

# Output dirs
BASE = Path('.')
RAW_DIR = BASE / 'proteins_raw'
RAW_DIR.mkdir(exist_ok=True, parents=True)
META_CSV = BASE / 'virus_metadata.csv'
TAX_CSV = BASE / 'virus_taxonomy.csv'
DOMAINS_DIR = BASE / 'domains'
DOMAINS_DIR.mkdir(exist_ok=True, parents=True)
TREES_DIR = BASE / 'trees'
TREES_DIR.mkdir(exist_ok=True, parents=True)


import time

def ncbi_esearch(db: str, term: str, retmax: int = 5) -> List[str]:
    params = {'db': db, 'term': term, 'retmode': 'json', 'retmax': retmax}
    for attempt in range(5):
        r = requests.get(f'{EUTILS}/esearch.fcgi', params=params, headers=UA, timeout=60)
        if r.status_code == 429:
            time.sleep(1.5 * (attempt + 1))
            continue
        r.raise_for_status()
        js = r.json()
        return js.get('esearchresult', {}).get('idlist', [])
    return []


def ncbi_efetch_fasta(db: str, ids: List[str]) -> Optional[str]:
    if not ids:
        return None
    params = {'db': db, 'id': ','.join(ids), 'rettype': 'fasta', 'retmode': 'text'}
    for attempt in range(5):
        r = requests.get(f'{EUTILS}/efetch.fcgi', params=params, headers=UA, timeout=60)
        if r.status_code == 429:
            time.sleep(1.5 * (attempt + 1))
            continue
        if r.status_code == 200 and r.text.strip():
            return r.text
        break
    return None


def ncbi_taxonomy_by_name(name: str) -> Tuple[Optional[str], Dict[str, str]]:
    # Try main name and aliases
    cand_names = [name] + ALIASES.get(name, [])
    taxid = None
    lineage = {}
    for nm in cand_names:
        ids = ncbi_esearch('taxonomy', f'"{nm}"[Scientific Name]')
        if not ids and 'virus' not in nm.lower():
            ids = ncbi_esearch('taxonomy', f'"{nm} virus"[All Names]')
        if ids:
            taxid = ids[0]
            break
    if not taxid:
        return None, {}
    # Fetch taxonomy XML and parse simple lineage with regex
    params = {'db': 'taxonomy', 'id': taxid}
    for attempt in range(5):
        r = requests.get(f'{EUTILS}/efetch.fcgi', params=params, headers=UA, timeout=60)
        if r.status_code == 429:
            time.sleep(1.5 * (attempt + 1))
            continue
        if r.status_code == 200:
            txt = r.text
            # extract lineage ranks using simple regex
            for rank in ['realm','kingdom','phylum','subphylum','class','order','family','genus','species']:
                m = re.search(fr'<Rank>{rank}</Rank>\s*<ScientificName>([^<]+)</ScientificName>', txt)
                if m:
                    lineage[rank] = m.group(1)
            m = re.search(r'<Lineage>([^<]+)</Lineage>', txt)
            if m:
                lineage['lineage'] = m.group(1)
            break
    return taxid, lineage


def pick_protein_term(virus: str, target: str) -> List[str]:
    if target == 'PB1':
        return [f'"{virus}"[Organism] AND (PB1 OR polymerase catalytic subunit OR RNA-dependent RNA polymerase PB1)']
    if target == 'PA':
        return [f'"{virus}"[Organism] AND (PA OR polymerase acidic OR endonuclease)']
    # L protein searches
    terms = [
        f'"{virus}"[Organism] AND (L[Gene Name] OR L protein OR RNA-dependent RNA polymerase OR polymerase L)'
    ]
    return terms


def fetch_protein(virus: str, target: str) -> Optional[Tuple[str, SeqRecord, str]]:
    # returns (accession, record, source)
    for term in pick_protein_term(virus, target):
        ids = ncbi_esearch('protein', term, retmax=5)
        if not ids:
            # try alias names
            for alias in ALIASES.get(virus, []):
                ids = ncbi_esearch('protein', term.replace(f'"{virus}"', f'"{alias}"'), retmax=5)
                if ids:
                    break
        if ids:
            fas = ncbi_efetch_fasta('protein', [ids[0]])
            if fas:
                from io import StringIO
                recs = list(SeqIO.parse(StringIO(fas), 'fasta'))
                if recs:
                    rec = recs[0]
                    acc = rec.id.split('|')[0]
                    return acc, rec, 'NCBI'
    # fallback: UniProt
    try:
        upq = None
        if target == 'PB1':
            upq = f'"{virus}" AND (gene:PB1)'
        elif target == 'PA':
            upq = f'"{virus}" AND (gene:PA)'
        else:
            upq = f'"{virus}" AND (gene:L OR name:polymerase)'
        url = 'https://rest.uniprot.org/uniprotkb/search'
        r = requests.get(url, params={'query': upq, 'format': 'fasta', 'size': 1}, headers=UA, timeout=60)
        if r.status_code == 200 and r.text.strip():
            from io import StringIO
            recs = list(SeqIO.parse(StringIO(r.text), 'fasta'))
            if recs:
                rec = recs[0]
                acc = rec.id.split('|')[1] if '|' in rec.id else rec.id
                return acc, rec, 'UniProt'
    except Exception:
        pass
    return None


def save_fasta(rec: SeqRecord, path: Path):
    with path.open('w') as fh:
        SeqIO.write(rec, fh, 'fasta')


def build_metadata_and_download():
    meta_rows = []
    for virus in VIRUSES:
        taxid, lineage = ncbi_taxonomy_by_name(virus)
        for target in PROTEIN_TARGETS[virus]:
            fetched = fetch_protein(virus, target)
            if not fetched:
                print(f'WARN: failed to fetch {virus} {target}')
                continue
            acc, rec, source = fetched
            fn = RAW_DIR / f"{virus.replace(' ','_')}_{target}_{acc}.fasta"
            save_fasta(rec, fn)
            row = {
                'Virus_name': virus,
                'Family': lineage.get('family',''),
                'Order': lineage.get('order',''),
                'Class': lineage.get('class',''),
                'Phylum': lineage.get('phylum',''),
                'NCBI_taxid': taxid or '',
                'Protein_name': target,
                'Accession_ID': acc,
                'Protein_length': len(rec.seq),
                'Downloaded_from': source,
            }
            meta_rows.append(row)
    # write CSV
    with META_CSV.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['Virus_name','Family','Order','Class','Phylum','NCBI_taxid','Protein_name','Accession_ID','Protein_length','Downloaded_from'])
        w.writeheader()
        for r in meta_rows:
            w.writerow(r)
    print(f'Metadata written: {META_CSV} with {len(meta_rows)} rows')


def find_endonuclease_domain(seq: str, virus: str, protein: str) -> Optional[Tuple[int,int]]:
    # Prefer N-terminus for cap-snatching endonucleases (Orthomyxo, Bunya, Arena)
    # Heuristic: if H-PD-D/E-K pattern occurs within first 250 aa, take 1-240; else take 1-240 for orthomyxo/bunya/arena; otherwise None
    upper = seq.upper()
    n = min(250, len(seq))
    window = upper[:n]
    if re.search(r'H.{0,20}P[D,E].{0,3}[D,E].{0,60}K', window):
        return (1, min(240, len(seq)))
    # decide based on taxonomic expectation via name
    if any(x in virus for x in ['Influenza','Thogoto','Dhori','Quaranfil','Hantaan','Seoul','Andes','Sin Nombre','Puumala','Dobrava','Tula','Prospect Hill','Rift Valley','La Crosse','Bunyamwera','Oropouche','Schmallenberg','Akabane','SFTS','Heartland','Uukuniemi','Toscana','Crimean-Congo','Nairobi sheep','Dugbe','Hazara','Lassa','LCMV','Junin','Machupo','Guanarito']):
        return (1, min(240, len(seq)))
    return None


def find_polymerase_palm(seq: str) -> Optional[Tuple[int,int]]:
    # Locate GDD motif and take +/-120 aa window
    upper = seq.upper()
    m = re.search(r'GDD', upper)
    if not m:
        m = re.search(r'GDNQ', upper)
    if not m:
        return None
    c = m.start()
    start = max(0, c-120)
    end = min(len(seq), c+80)
    if end-start < 80:
        return None
    return (start+1, end)


def load_all_sequences() -> List[Tuple[str,str,str,str]]:
    # returns list of (label, virus, protein, sequence)
    items = []
    for fa in RAW_DIR.glob('*.fasta'):
        rec = next(SeqIO.parse(str(fa), 'fasta'))
        parts = fa.stem.split('_')
        protein = parts[-2]
        virus = ' '.join(parts[:-2]).replace('_',' ')
        label = f"{virus}|{protein}"
        items.append((label, virus, protein, str(rec.seq)))
    return items


def make_domain_fastas():
    items = load_all_sequences()
    # Endonuclease set
    endo_records = []
    for label, virus, prot, seq in items:
        if prot not in ('PA','L'):
            continue
        dom = find_endonuclease_domain(seq, virus, prot)
        if dom is None:
            continue
        s,e = dom
        sub = seq[s-1:e]
        rec = SeqRecord(Seq(sub), id=label.replace(' ','_'), description=f'endonuclease {s}-{e}')
        endo_records.append(rec)
    if endo_records:
        endo_fa = DOMAINS_DIR / 'endonuclease_domains.fasta'
        SeqIO.write(endo_records, str(endo_fa), 'fasta')
        print('Wrote', endo_fa, len(endo_records))
    # Polymerase palm set
    pol_records = []
    for label, virus, prot, seq in items:
        if prot not in ('PB1','L'):
            continue
        dom = find_polymerase_palm(seq)
        if dom is None:
            continue
        s,e = dom
        sub = seq[s-1:e]
        rec = SeqRecord(Seq(sub), id=label.replace(' ','_'), description=f'palm {s}-{e}')
        pol_records.append(rec)
    if pol_records:
        pol_fa = DOMAINS_DIR / 'polymerase_palm_domains.fasta'
        SeqIO.write(pol_records, str(pol_fa), 'fasta')
        print('Wrote', pol_fa, len(pol_records))


def align_and_trim(in_fa: Path, out_prefix: str):
    aligned = DOMAINS_DIR / f'{out_prefix}_aligned.fasta'
    trimmed = DOMAINS_DIR / f'{out_prefix}_trimmed.fasta'
    os.system(f"mafft --auto --thread -1 {in_fa} > {aligned} 2> {DOMAINS_DIR/(out_prefix + '_mafft.log')}")
    # simple trimming: drop columns with >50% gaps
    aln = AlignIO.read(str(aligned), 'fasta')
    keep_cols = []
    for i in range(aln.get_alignment_length()):
        col = aln[:, i]
        gaps = col.count('-') + col.count('.')
        if gaps/len(aln) <= 0.5:
            keep_cols.append(i)
    trimmed_aln = None
    for idx in keep_cols:
        if trimmed_aln is None:
            trimmed_aln = aln[:, idx:idx+1]
        else:
            trimmed_aln += aln[:, idx:idx+1]
    if trimmed_aln is None:
        trimmed_aln = aln
    AlignIO.write(trimmed_aln, str(trimmed), 'fasta')
    # NJ tree
    calc = DistanceCalculator('blosum62')
    dm = calc.get_distance(trimmed_aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    newick = TREES_DIR / f'{out_prefix}.nwk'
    Phylo.write(tree, str(newick), 'newick')
    # figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10,7), dpi=150)
    ax = fig.add_subplot(1,1,1)
    Phylo.draw(tree, do_show=False, axes=ax)
    fig.tight_layout()
    fig.savefig(TREES_DIR / f'{out_prefix}.png')
    print('Tree written', newick)


def main():
    build_metadata_and_download()
    make_domain_fastas()
    # Align and tree building for two domain sets
    endo_fa = DOMAINS_DIR / 'endonuclease_domains.fasta'
    if endo_fa.exists():
        align_and_trim(endo_fa, 'endonuclease_domains')
    pol_fa = DOMAINS_DIR / 'polymerase_palm_domains.fasta'
    if pol_fa.exists():
        align_and_trim(pol_fa, 'polymerase_palm_domains')

if __name__ == '__main__':
    main()
