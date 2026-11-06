#!/usr/bin/env python3
import os
import sys
import time
import json
import math
import csv
from typing import List, Dict, Tuple, Set

import requests
import pandas as pd
from io import StringIO

# BioPython for NCBI parsing
from Bio import Entrez, SeqIO

# Configuration
USER_AGENT = "PlayMoleculeAI/1.0 (Acellera; data-retrieval; contact: support@acellera.com)"
Entrez.email = os.environ.get("NCBI_EMAIL", "playmoleculeai@acellera.com")
Entrez.api_key = os.environ.get("NCBI_API_KEY")

BASE_OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "data"))
UNIPROT_DIR = os.path.join(BASE_OUT, "uniprot")
NCBI_DIR = os.path.join(BASE_OUT, "ncbi")
PDB_DIR = os.path.join(BASE_OUT, "pdb")
for d in (UNIPROT_DIR, NCBI_DIR, PDB_DIR):
    os.makedirs(d, exist_ok=True)

# Helpers

def get(url, **kwargs):
    headers = kwargs.pop("headers", {})
    headers.setdefault("User-Agent", USER_AGENT)
    r = requests.get(url, headers=headers, timeout=120, **kwargs)
    if r.status_code == 429:
        # simple backoff
        time.sleep(2)
        r = requests.get(url, headers=headers, timeout=120, **kwargs)
    r.raise_for_status()
    return r


def uniprot_download_all():
    # Query both current species name and common synonym
    query = '((organism_name:"Andes orthohantavirus") OR (organism_name:"Andes virus"))'

    # Metadata TSV
    fields = [
        "accession",
        "id",
        "protein_name",
        "gene_primary",
        "gene_names",
        "organism_name",
        "organism_id",
        "length",
        "reviewed",
        "cc_subcellular_location",
        "cc_function",
        "xref_pdb",
        "protein_existence",
    ]
    url_tsv = (
        "https://rest.uniprot.org/uniprotkb/stream?" +
        f"query={requests.utils.quote(query)}&format=tsv&fields={','.join(fields)}"
    )
    tsv_path = os.path.join(UNIPROT_DIR, "andes_orthohantavirus_uniprot.tsv")
    print(f"Downloading UniProt metadata TSV -> {tsv_path}")
    r = get(url_tsv)
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write(r.text)

    # Load into DataFrame and save CSV cleanly
    df = pd.read_csv(tsv_path, sep="\t")
    # Normalize PDB IDs column to a pipe-separated list without spaces
    pdb_source_col = None
    for cand in ["Cross-reference (PDB)", "PDB"]:
        if cand in df.columns:
            pdb_source_col = cand
            break
    if pdb_source_col:
        df["pdb_ids"] = df[pdb_source_col].fillna("").astype(str).str.replace("; ", ";", regex=False)
        df["pdb_ids"] = df["pdb_ids"].str.replace(" ", "", regex=False)
    else:
        df["pdb_ids"] = ""

    csv_path = os.path.join(UNIPROT_DIR, "andes_orthohantavirus_uniprot_metadata.csv")
    df.to_csv(csv_path, index=False)

    # FASTA
    url_fasta = (
        "https://rest.uniprot.org/uniprotkb/stream?" +
        f"query={requests.utils.quote(query)}&format=fasta"
    )
    fasta_path = os.path.join(UNIPROT_DIR, "andes_orthohantavirus_uniprot.fasta")
    print(f"Downloading UniProt FASTA -> {fasta_path}")
    r = get(url_fasta)
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(r.text)

    # Return set of PDB IDs
    pdb_ids: Set[str] = set()
    if "pdb_ids" in df.columns:
        for ids in df["pdb_ids"].dropna().astype(str):
            for pid in str(ids).replace(";", ",").split(","):
                pid = pid.strip()
                if pid:
                    # UniProt may include chain or other annotations: e.g., 3XYZ
                    # Filter to 4-char PDB IDs
                    if len(pid) == 4 and pid.isalnum():
                        pdb_ids.add(pid.upper())
    return pdb_ids, csv_path, fasta_path


def ncbi_find_taxids() -> List[str]:
    # Search NCBI Taxonomy for exact scientific name matches
    term = '"Andes orthohantavirus"[SCIN] OR "Andes virus"[SCIN]'
    h = Entrez.esearch(db="taxonomy", term=term)  # XML by default
    res = Entrez.read(h)
    taxids = res.get("IdList", [])
    print(f"Found taxonomy IDs: {taxids}")
    return taxids


def ncbi_protein_ids_by_organism_name() -> List[str]:
    term = '"Andes orthohantavirus"[Organism] OR "Andes virus"[Organism]'
    ids: List[str] = []
    retstart = 0
    retmax = 5000
    while True:
        h = Entrez.esearch(db="protein", term=term, retstart=retstart, retmax=retmax)
        res = Entrez.read(h)
        batch = res.get("IdList", [])
        if not batch:
            break
        ids.extend(batch)
        retstart += len(batch)
        count = int(res.get("Count", "0"))
        print(f"NCBI protein IDs (name search): {len(ids)}/{count}")
        if len(ids) >= count:
            break
        time.sleep(0.34)
    return ids


def ncbi_protein_ids_for_taxid(taxid: str) -> List[str]:
    # Use organism expansion to include all children (strains/isolates)
    term = f"txid{taxid}[Organism:exp]"
    ids: List[str] = []
    retstart = 0
    retmax = 5000
    while True:
        h = Entrez.esearch(db="protein", term=term, retstart=retstart, retmax=retmax)
        res = Entrez.read(h)
        batch = res.get("IdList", [])
        if not batch:
            break
        ids.extend(batch)
        retstart += len(batch)
        count = int(res.get("Count", "0"))
        print(f"NCBI protein IDs: {len(ids)}/{count}")
        if len(ids) >= count:
            break
        time.sleep(0.34)  # polite delay (~3 req/s)
    return ids


def ncbi_fetch_fasta(ids: List[str], out_fasta: str):
    print(f"Fetching NCBI protein FASTA for {len(ids)} IDs -> {out_fasta}")
    with open(out_fasta, "w", encoding="utf-8") as out:
        for i in range(0, len(ids), 500):
            batch = ids[i:i+500]
            h = Entrez.efetch(db="protein", id=",".join(batch), rettype="fasta", retmode="text")
            data = h.read()
            out.write(data)
            time.sleep(0.34)


def ncbi_fetch_metadata(ids: List[str], out_csv: str):
    print(f"Fetching NCBI protein metadata for {len(ids)} IDs -> {out_csv}")
    records = []
    for i in range(0, len(ids), 200):
        batch = ids[i:i+200]
        h = Entrez.efetch(db="protein", id=",".join(batch), rettype="gb", retmode="text")
        for rec in SeqIO.parse(h, "genbank"):
            # Extract key metadata
            acc = rec.id
            descr = rec.description
            length = len(rec.seq)
            organism = rec.annotations.get("organism")
            # Try to get strain/isolate from source features
            strain = isolate = host = collection_date = country = None
            for feat in rec.features:
                if feat.type == "source":
                    strain = feat.qualifiers.get("strain", [None])[0]
                    isolate = feat.qualifiers.get("isolate", [None])[0]
                    host = feat.qualifiers.get("host", [None])[0]
                    collection_date = feat.qualifiers.get("collection_date", [None])[0]
                    country = feat.qualifiers.get("country", [None])[0]
            product = gene = None
            for feat in rec.features:
                if feat.type == "CDS":
                    if not product:
                        product = feat.qualifiers.get("product", [None])[0]
                    if not gene:
                        gene = feat.qualifiers.get("gene", [None])[0]
            records.append({
                "accession": acc,
                "description": descr,
                "length": length,
                "organism": organism,
                "strain": strain,
                "isolate": isolate,
                "host": host,
                "collection_date": collection_date,
                "country": country,
                "gene": gene,
                "product": product,
            })
        time.sleep(0.34)
    # Save CSV
    df = pd.DataFrame.from_records(records)
    df.to_csv(out_csv, index=False)


def rcsb_entry_metadata(pdb_id: str) -> Dict:
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        r = get(url)
        return r.json()
    except Exception as e:
        print(f"Failed metadata for {pdb_id}: {e}")
        return {}


def download_pdb(pdb_id: str) -> Tuple[str, Dict]:
    # Try biological assembly PDB first, then asymmetric unit PDB
    out_pdb1 = os.path.join(PDB_DIR, f"{pdb_id}_assembly1.pdb")
    out_pdb = os.path.join(PDB_DIR, f"{pdb_id}.pdb")
    ok = False
    try:
        url1 = f"https://files.rcsb.org/download/{pdb_id}.pdb1"
        r = requests.get(url1, headers={"User-Agent": USER_AGENT}, timeout=120)
        if r.status_code == 200 and r.text.strip():
            with open(out_pdb1, "w", encoding="utf-8") as f:
                f.write(r.text)
            ok = True
    except Exception:
        pass
    if not ok:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        try:
            r = get(url)
            with open(out_pdb, "w", encoding="utf-8") as f:
                f.write(r.text)
        except requests.HTTPError as e:
            # Fallback to mmCIF if legacy PDB not available
            cif_path = os.path.join(PDB_DIR, f"{pdb_id}.cif")
            url_cif = f"https://files.rcsb.org/download/{pdb_id}.cif"
            rc = get(url_cif)
            with open(cif_path, "w", encoding="utf-8") as f:
                f.write(rc.text)
            out_pdb = cif_path
    meta = rcsb_entry_metadata(pdb_id)
    return (out_pdb1 if ok else out_pdb), meta


def save_pdb_metadata(all_meta: Dict[str, Dict], out_csv: str):
    rows = []
    for pid, m in all_meta.items():
        title = m.get("struct", {}).get("title") if m else None
        exptl = None
        if m and m.get("exptl"):
            methods = [e.get("method") for e in m.get("exptl") if isinstance(e, dict)]
            exptl = ";".join([x for x in methods if x]) if methods else None
        res = None
        if m:
            info = m.get("rcsb_entry_info", {})
            rc = info.get("resolution_combined")
            if isinstance(rc, list) and rc:
                res = rc[0]
            elif isinstance(rc, (int, float)):
                res = rc
        dep = m.get("rcsb_accession_info", {}).get("deposit_date") if m else None
        rel = m.get("rcsb_accession_info", {}).get("initial_release_date") if m else None
        polymer_entities = m.get("polymer_entity_count") if m else None
        nonpolymer_entities = m.get("nonpolymer_entity_count") if m else None
        assemblies = m.get("rcsb_entry_container_identifiers", {}).get("assembly_ids") if m else None
        rows.append({
            "pdb_id": pid,
            "title": title,
            "experimental_method": exptl,
            "resolution_A": res,
            "deposit_date": dep,
            "release_date": rel,
            "polymer_entity_count": polymer_entities,
            "nonpolymer_entity_count": nonpolymer_entities,
            "assembly_ids": ";".join(assemblies) if assemblies else None,
        })
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def main():
    print("== UniProt retrieval ==")
    pdb_ids_uniprot, uniprot_csv, uniprot_fasta = uniprot_download_all()

    print("== NCBI retrieval ==")
    taxids = ncbi_find_taxids()
    all_ncbi_ids: Set[str] = set()
    for tx in taxids:
        ids = ncbi_protein_ids_for_taxid(tx)
        all_ncbi_ids.update(ids)
    # Fallback to organism-name search if taxonomy IDs not found
    if not all_ncbi_ids:
        print("No taxonomy-based NCBI IDs found; falling back to organism name search.")
        ids_by_name = ncbi_protein_ids_by_organism_name()
        all_ncbi_ids.update(ids_by_name)
    all_ncbi_ids = list(all_ncbi_ids)
    print(f"Total unique NCBI protein IDs: {len(all_ncbi_ids)}")

    ncbi_fasta = os.path.join(NCBI_DIR, "andes_orthohantavirus_ncbi.fasta")
    ncbi_meta_csv = os.path.join(NCBI_DIR, "andes_orthohantavirus_ncbi_metadata.csv")
    if all_ncbi_ids:
        ncbi_fetch_fasta(all_ncbi_ids, ncbi_fasta)
        ncbi_fetch_metadata(all_ncbi_ids, ncbi_meta_csv)

    print("== PDB download from UniProt cross-references ==")
    all_pdb_meta: Dict[str, Dict] = {}
    downloaded = []
    for i, pid in enumerate(sorted(pdb_ids_uniprot)):
        print(f"[{i+1}/{len(pdb_ids_uniprot)}] Downloading {pid}")
        path, meta = download_pdb(pid)
        downloaded.append((pid, path))
        all_pdb_meta[pid] = meta
        time.sleep(0.2)

    pdb_meta_csv = os.path.join(PDB_DIR, "pdb_metadata.csv")
    save_pdb_metadata(all_pdb_meta, pdb_meta_csv)

    # Manifest JSON with paths
    manifest = {
        "uniprot": {
            "metadata_csv": os.path.relpath(uniprot_csv, BASE_OUT),
            "fasta": os.path.relpath(uniprot_fasta, BASE_OUT),
            "count_entries": int(pd.read_csv(uniprot_csv).shape[0]),
        },
        "ncbi": {
            "metadata_csv": os.path.relpath(ncbi_meta_csv, BASE_OUT) if all_ncbi_ids else None,
            "fasta": os.path.relpath(ncbi_fasta, BASE_OUT) if all_ncbi_ids else None,
            "count_entries": len(all_ncbi_ids),
        },
        "pdb": {
            "downloaded": [
                {"pdb_id": pid, "path": os.path.relpath(pth, BASE_OUT)} for pid, pth in downloaded
            ],
            "metadata_csv": os.path.relpath(pdb_meta_csv, BASE_OUT) if downloaded else None,
        }
    }
    manifest_path = os.path.join(BASE_OUT, "manifest_andes_orthohantavirus.json")
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print("\nDone. Outputs in:")
    print(f" - {UNIPROT_DIR}")
    print(f" - {NCBI_DIR}")
    print(f" - {PDB_DIR}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
