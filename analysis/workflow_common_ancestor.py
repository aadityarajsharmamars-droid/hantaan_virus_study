import os
import re
import json
import csv
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import requests
from Bio import SeqIO, AlignIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import Consensus

UA={'User-Agent':'PlayMoleculeAI/Negarnaviricota (contact:user)'}
EUTILS='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

VIRUSES=[
"Andes virus","Hantaan virus","Seoul virus","Sin Nombre virus","Puumala virus","Dobrava-Belgrade virus","Tula virus","Prospect Hill virus",
"La Crosse virus","California encephalitis virus","Bunyamwera virus","Oropouche virus","Schmallenberg virus","Akabane virus","Rift Valley fever virus","SFTS virus","Heartland virus","Uukuniemi virus","Toscana virus","Crimean-Congo hemorrhagic fever virus","Nairobi sheep disease virus","Dugbe virus","Hazara virus",
"Lassa virus","LCMV","Junin virus","Machupo virus","Guanarito virus",
"Influenza A virus","Influenza B virus","Influenza C virus","Influenza D virus","Thogoto virus","Dhori virus","Quaranfil virus","Infectious salmon anemia virus",
"Vesicular stomatitis virus","Rabies virus","Measles virus","Ebola virus"
]

# Map families/orders for protein choice guidance
FAMILY_HINTS={
    'Hantaan virus':'Hantaviridae', 'Andes virus':'Hantaviridae','Seoul virus':'Hantaviridae','Sin Nombre virus':'Hantaviridae','Puumala virus':'Hantaviridae','Dobrava-Belgrade virus':'Hantaviridae','Tula virus':'Hantaviridae','Prospect Hill virus':'Hantaviridae',
    'La Crosse virus':'Peribunyaviridae','California encephalitis virus':'Peribunyaviridae','Bunyamwera virus':'Peribunyaviridae','Oropouche virus':'Peribunyaviridae','Schmallenberg virus':'Peribunyaviridae','Akabane virus':'Peribunyaviridae',
    'Rift Valley fever virus':'Phenuiviridae','SFTS virus':'Phenuiviridae','Heartland virus':'Phenuiviridae','Uukuniemi virus':'Phenuiviridae','Toscana virus':'Phenuiviridae',
    'Crimean-Congo hemorrhagic fever virus':'Nairoviridae','Nairobi sheep disease virus':'Nairoviridae','Dugbe virus':'Nairoviridae','Hazara virus':'Nairoviridae',
    'Influenza A virus':'Orthomyxoviridae','Influenza B virus':'Orthomyxoviridae','Influenza C virus':'Orthomyxoviridae','Influenza D virus':'Orthomyxoviridae','Thogoto virus':'Orthomyxoviridae','Dhori virus':'Orthomyxoviridae','Quaranfil virus':'Orthomyxoviridae','Infectious salmon anemia virus':'Orthomyxoviridae',
    'Lassa virus':'Arenaviridae','LCMV':'Arenaviridae','Junin virus':'Arenaviridae','Machupo virus':'Arenaviridae','Guanarito virus':'Arenaviridae',
    'Vesicular stomatitis virus':'Rhabdoviridae','Rabies virus':'Rhabdoviridae','Measles virus':'Paramyxoviridae','Ebola virus':'Filoviridae'
}

# Which proteins to fetch per family/order rule
PROTEIN_RULES={
    'Orthomyxoviridae':['PA','PB1'],
    'default_L':['L']
}

project=Path('.')
raw_dir=project/'proteins_raw'
raw_dir.mkdir(parents=True, exist_ok=True)

meta_rows=[]


def ncbi_taxid(name:str)->Optional[str]:
    r=requests.get(EUTILS+'esearch.fcgi', params={'db':'taxonomy','term':name,'retmode':'json'}, headers=UA, timeout=60)
    if r.status_code==200:
        ids=r.json()['esearchresult']['idlist']
        return ids[0] if ids else None
    return None

def ncbi_taxonomy_lineage(taxid:str)->Dict[str,str]:
    import xml.etree.ElementTree as ET
    r=requests.get(EUTILS+'efetch.fcgi', params={'db':'taxonomy','id':taxid,'retmode':'xml'}, headers=UA, timeout=60)
    text=r.text
    if not text.strip().startswith('<?xml') and not text.strip().startswith('<'):
        # fallback to esummary
        r2=requests.get(EUTILS+'esummary.fcgi', params={'db':'taxonomy','id':taxid,'retmode':'json'}, headers=UA, timeout=60)
        js=r2.json()
        doc=js['result'][taxid]
        lineage=doc.get('lineage','')
        # rough parse
        mapping={'Realm':'Riboviria' if 'Riboviria' in lineage else '', 'Kingdom':'', 'Phylum':'', 'Class':'', 'Order':'', 'Family':'', 'Genus':'', 'Species':doc.get('scientificname','')}
        return mapping
    root=ET.fromstring(text)
    taxon=root.find('.//Taxon')
    ranks={}
    lineage_list=[]
    for node in taxon.findall('LineageEx/Taxon'):  # nested Taxon elements
        rank=(node.findtext('Rank') or '').lower()
        name=node.findtext('ScientificName') or ''
        if rank:
            ranks[rank]=name
            lineage_list.append((rank,name))
    sci_name=taxon.findtext('ScientificName') or ''
    lineage_str=taxon.findtext('Lineage') or ''
    mapping={
        'Realm': 'Riboviria' if 'Riboviria' in lineage_str or sci_name=='Riboviria' else '',
        'Kingdom': '',
        'Phylum': ranks.get('phylum',''),
        'Class': ranks.get('class',''),
        'Order': ranks.get('order',''),
        'Family': ranks.get('family',''),
        'Genus': ranks.get('genus',''),
        'Species': sci_name
    }
    return mapping


def esearch_protein(query:str, retmax:int=20)->List[str]:
    r=requests.get(EUTILS+'esearch.fcgi', params={'db':'protein','term':query,'retmode':'json','retmax':retmax}, headers=UA, timeout=60)
    if r.status_code!=200:
        return []
    return r.json()['esearchresult'].get('idlist',[])


def efetch_protein_fasta(acc:str)->Optional[str]:
    r=requests.get(EUTILS+'efetch.fcgi', params={'db':'protein','id':acc,'rettype':'fasta','retmode':'text'}, headers=UA, timeout=60)
    if r.status_code==200 and r.text and r.text.strip().startswith('>'):
        return r.text
    return None

# UniProt fallback
UNIPROT_BASE='https://rest.uniprot.org/uniprotkb/search'

def uniprot_fetch_one(virus:str, protein_keyword:str)->Optional[Tuple[str,str]]:
    # returns (fasta, accession)
    query=f'(organism_name:"{virus}") AND ({protein_keyword})'
    for q in [query+" AND reviewed:true", query]:
        r=requests.get(UNIPROT_BASE, params={'query':q,'format':'fasta','size':1}, headers=UA, timeout=60)
        if r.status_code==200 and r.text.strip().startswith('>'):
            header=r.text.splitlines()[0]
            m=re.search(r'>\w+\|(\w+)\|', header)
            acc=m.group(1) if m else 'UNIPROT'
            return r.text, acc
    return None


def esummary_protein(acc:str)->Dict:
    r=requests.get(EUTILS+'esummary.fcgi', params={'db':'protein','id':acc,'retmode':'json'}, headers=UA, timeout=60)
    if r.status_code==200:
        return r.json()['result'].get(acc, {})
    return {}


def pick_proteins(virus:str, fam_hint:str)->List[Tuple[str,str,str,str]]:
    # returns list of (protein_name, accession, source, filepath)
    picks=[]
    vtag=virus.replace(' ','_')
    if fam_hint=='Orthomyxoviridae':
        for gene in ['PA','PB1']:
            # Try NCBI first
            q=f'({virus}[Organism]) AND ({gene}[Gene] OR {gene}[Protein Name]) NOT partial'
            ids=esearch_protein(q, retmax=5)
            if ids:
                acc=ids[0]
                fasta=efetch_protein_fasta(acc)
                if fasta:
                    path=raw_dir/f"{vtag}__{gene}__{acc}.fasta"
                    path.write_text(fasta)
                    picks.append((gene, acc, 'NCBI Protein', str(path)))
                    continue
            # Fallback UniProt
            key = f'gene_exact:{gene} OR protein_name:"RNA-directed RNA polymerase catalytic subunit"' if gene=='PB1' else 'protein_name:endonuclease OR gene_exact:PA'
            up=uniprot_fetch_one(virus, key)
            if up:
                fasta, uacc=up
                path=raw_dir/f"{vtag}__{gene}__{uacc}.fasta"
                path.write_text(fasta)
                picks.append((gene, uacc, 'UniProt', str(path)))
    else:
        # L protein
        q=f'({virus}[Organism]) AND (L[Gene] OR polymerase L[Protein Name] OR RNA-dependent RNA polymerase[Protein Name]) NOT partial'
        ids=esearch_protein(q, retmax=5)
        if ids:
            acc=ids[0]
            fasta=efetch_protein_fasta(acc)
            if fasta:
                path=raw_dir/f"{vtag}__L__{acc}.fasta"
                path.write_text(fasta)
                picks.append(('L', acc, 'NCBI Protein', str(path)))
        if not picks:
            up=uniprot_fetch_one(virus, 'protein_name:(polymerase L) OR gene_exact:L OR protein_name:"RNA-directed RNA polymerase"')
            if up:
                fasta,uacc=up
                path=raw_dir/f"{vtag}__L__{uacc}.fasta"
                path.write_text(fasta)
                picks.append(('L', uacc, 'UniProt', str(path)))
    return picks


def ensure_fasta_written(virus:str, prot_name:str, acc:str)->Tuple[str,int]:
    fasta=efetch_protein_fasta(acc)
    if not fasta:
        raise RuntimeError(f'No FASTA for {virus} {prot_name} {acc}')
    # parse length
    from io import StringIO
    seq=None
    handle=StringIO(fasta)
    for rec in SeqIO.parse(handle, 'fasta'):
        seq=str(rec.seq)
        break
    if seq is None:
        # fallback parse manually
        lines=fasta.strip().split('\n')
        seq=''.join(l.strip() for l in lines if not l.startswith('>'))
    length=len(seq)
    out_path=raw_dir/f"{virus.replace(' ','_')}__{prot_name}__{acc}.fasta"
    out_path.write_text(fasta)
    return str(out_path), length


def build_metadata():
    for virus in VIRUSES:
        taxid=ncbi_taxid(virus)
        lineage=ncbi_taxonomy_lineage(taxid) if taxid else {k:'' for k in ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species']}
        fam_hint=FAMILY_HINTS.get(virus, '')
        picks=pick_proteins(virus, fam_hint)
        if not picks:
            print('WARN no proteins for', virus)
            continue
        for prot_name, acc, source, path in picks:
            try:
                # compute length by parsing
                from io import StringIO
                txt=Path(path).read_text()
                handle=StringIO(txt)
                rec=next(SeqIO.parse(handle,'fasta'))
                length=len(rec.seq)
            except Exception as e:
                print('ERR parsing', virus, prot_name, acc, e)
                continue
            meta_rows.append({
                'Virus_name': virus,
                'Family': lineage.get('Family',''),
                'Order': lineage.get('Order',''),
                'Class': lineage.get('Class',''),
                'Phylum': lineage.get('Phylum',''),
                'Realm': lineage.get('Realm',''),
                'NCBI_taxid': taxid or '',
                'Protein_name': prot_name,
                'Accession_ID': acc,
                'Protein_length': length,
                'Downloaded_from': source,
                'FASTA_path': path
            })
    # write metadata
    with open(project/'virus_metadata.csv','w', newline='') as fh:
        w=csv.DictWriter(fh, fieldnames=['Virus_name','Family','Order','Class','Phylum','Realm','NCBI_taxid','Protein_name','Accession_ID','Protein_length','Downloaded_from','FASTA_path'])
        w.writeheader()
        for r in meta_rows:
            w.writerow(r)


def parse_fasta_file(path:Path)->SeqRecord:
    return next(SeqIO.parse(str(path), 'fasta'))


def extract_endonuclease_sequences()->Path:
    # build endonuclease_domains.fasta from N-termini for cap-snatching families (Orthomyxoviridae, Bunyavirales*, Arenaviridae)
    families_cap={'Orthomyxoviridae','Hantaviridae','Peribunyaviridae','Phenuiviridae','Nairoviridae','Tospoviridae','Arenaviridae'}
    records=[]
    for r in meta_rows:
        if r['Family'] not in families_cap:
            continue
        if r['Protein_name'] in ('PA','L'):
            rec=parse_fasta_file(Path(r['FASTA_path']))
            seq=str(rec.seq).replace('U','X')
            dom=seq[:250] if len(seq)>=150 else seq
            records.append(SeqRecord(Seq(dom), id=f"{r['Virus_name'].replace(' ','_')}__{r['Protein_name']}", description=''))
    out=project/'endonuclease_domains.fasta'
    SeqIO.write(records, str(out), 'fasta')
    return out


def extract_polymerase_core()->Path:
    # search for GDD motif and extract +/-150 aa window
    records=[]
    for r in meta_rows:
        if r['Protein_name'] in ('PB1','L'):
            rec=parse_fasta_file(Path(r['FASTA_path']))
            seq=str(rec.seq)
            m=re.search(r'[ASG]DD', seq)
            if m:
                i=m.start()
                s=max(0, i-150)
                e=min(len(seq), i+150)
                core=seq[s:e]
            else:
                # fallback C-terminal 300aa
                core=seq[-300:]
            records.append(SeqRecord(Seq(core), id=f"{r['Virus_name'].replace(' ','_')}__{r['Protein_name']}", description=''))
    out=project/'polymerase_core.fasta'
    SeqIO.write(records, str(out), 'fasta')
    return out


def run_mafft(in_fa:Path, out_fa:Path):
    os.system(f"mafft --maxiterate 1000 --localpair {in_fa} > {out_fa} 2> mafft_{in_fa.stem}.log")


def trim_alignment(in_fa:Path, out_fa:Path, gap_thresh:float=0.5):
    aln=AlignIO.read(str(in_fa), 'fasta')
    keep_cols=[]
    for i in range(aln.get_alignment_length()):
        col=aln[:, i]
        gaps=col.count('-') + col.count('.')
        if gaps/len(aln) <= gap_thresh:
            keep_cols.append(i)
    trimmed=None
    for idx in keep_cols:
        trimmed = aln[:, idx:idx+1] if trimmed is None else trimmed + aln[:, idx:idx+1]
    AlignIO.write(trimmed, str(out_fa), 'fasta')


def build_nj_tree(aln_fa:Path, outprefix:str, boot:int=100):
    aln=AlignIO.read(str(aln_fa), 'fasta')
    calc=DistanceCalculator('blosum62')
    constructor=DistanceTreeConstructor(calc, method='nj')
    base_tree=constructor.build_tree(aln)
    Phylo.write(base_tree, f"{outprefix}.nwk", 'newick')
    # Bootstrap consensus
    cons=Consensus.bootstrap_consensus(aln, boot, constructor, Consensus.majority_consensus)
    Phylo.write(cons, f"{outprefix}.bootstrap.nwk", 'newick')
    # Plot
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(10,8), dpi=160)
    ax=fig.add_subplot(1,1,1)
    Phylo.draw(cons, do_show=False, axes=ax)
    fig.tight_layout()
    fig.savefig(f"{outprefix}.png")


def build_taxonomy_tree()->Path:
    # Create a simple taxonomy tree by stitching shared nodes from lineages
    # We'll build a nested dict root->...->species and then convert to Newick
    # Fetch full lineage for each species again to ensure completeness
    def lineage_list(taxid):
        import xml.etree.ElementTree as ET
        r=requests.get(EUTILS+'efetch.fcgi', params={'db':'taxonomy','id':taxid,'retmode':'xml'}, headers=UA, timeout=60)
        txt=r.text
        try:
            root=ET.fromstring(txt)
        except Exception:
            # fallback simple list with species only
            return [('species', taxid)]
        taxon=root.find('.//Taxon')
        lst=[]
        for node in taxon.findall('LineageEx/Taxon'):
            rank=(node.findtext('Rank') or '').lower()
            name=node.findtext('ScientificName') or ''
            if rank in ('superkingdom','phylum','class','order','family','genus'):
                lst.append((rank,name))
        lst.append(('species', taxon.findtext('ScientificName') or ''))
        return lst
    tree_dict={}
    for v in VIRUSES:
        taxid=ncbi_taxid(v)
        if not taxid: continue
        lineage= lineage_list(taxid)
        d=tree_dict
        for rank,name in lineage:
            d=d.setdefault(name,{})
    # convert dict to Newick
    def to_newick(d):
        if not d:
            return ''
        parts=[]
        for name,child in d.items():
            sub=to_newick(child)
            parts.append((name, sub))
        return '('+','.join([ (p[1]+p[0]) if p[1] else p[0] for p in parts ])+')'
    newick=to_newick(tree_dict)+';'
    out=project/'taxonomy_tree.nwk'
    Path(out).write_text(newick)
    return out


def main():
    build_metadata()
    # Endonuclease set
    endo_fa=extract_endonuclease_sequences()
    run_mafft(endo_fa, project/'aligned_endonuclease.fasta')
    trim_alignment(project/'aligned_endonuclease.fasta', project/'trimmed_endonuclease.fasta', 0.5)
    build_nj_tree(project/'trimmed_endonuclease.fasta', 'tree_endonuclease', boot=100)

    # Polymerase core set
    pol_fa=extract_polymerase_core()
    run_mafft(pol_fa, project/'aligned_polymerase.fasta')
    trim_alignment(project/'aligned_polymerase.fasta', project/'trimmed_polymerase.fasta', 0.5)
    build_nj_tree(project/'trimmed_polymerase.fasta', 'tree_polymerase', boot=100)

    # Taxonomy tree
    build_taxonomy_tree()

    print('Workflow complete.')

if __name__=='__main__':
    main()
