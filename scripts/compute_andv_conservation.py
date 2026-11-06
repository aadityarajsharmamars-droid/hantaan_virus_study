#!/usr/bin/env python3
import os
import re
import json
import math
import tempfile
import subprocess
from collections import Counter, defaultdict
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
DATADIR = os.path.join(BASE, "data")
CONSDIR = os.path.join(DATADIR, "conservation")
os.makedirs(CONSDIR, exist_ok=True)

NCBI_FASTA = os.path.join(DATADIR, "ncbi", "andes_orthohantavirus_ncbi.fasta")
NCBI_META = os.path.join(DATADIR, "ncbi", "andes_orthohantavirus_ncbi_metadata.csv")
UNIPROT_FASTA = os.path.join(DATADIR, "uniprot", "andes_orthohantavirus_uniprot.fasta")

# Reference UniProt accessions and lengths
REFS = {
    "N": {"acc": "O36307", "length": 428},
    "GPC": {"acc": "Q9E006", "length": 1138, "gn_range": (19, 651), "gc_range": (652, 1138)},
    "L": {"acc": "Q9E005", "length": 2153},
}

AA20 = list("ACDEFGHIKLMNPQRSWTYV")
AASET = set(AA20)
LOG2_20 = math.log2(20.0)


def read_fasta_dict(path: str) -> Dict[str, str]:
    seqs = {}
    for rec in SeqIO.parse(path, "fasta"):
        # accession is first token before whitespace
        header = rec.description
        acc = header.split()[0]
        if acc.startswith('>'):
            acc = acc[1:]
        seqs[acc] = str(rec.seq).upper()
    return seqs


def get_uniprot_seq(acc: str) -> Tuple[str, str]:
    """Return (header, sequence) for a UniProt accession from the UniProt FASTA file."""
    for rec in SeqIO.parse(UNIPROT_FASTA, "fasta"):
        if f"|{acc}|" in rec.description or rec.id.endswith(acc):
            return (rec.description, str(rec.seq).upper())
    raise RuntimeError(f"UniProt accession {acc} not found in FASTA {UNIPROT_FASTA}")


def classify_ncbi_records() -> pd.DataFrame:
    df = pd.read_csv(NCBI_META)
    # Normalize fields to lowercase strings
    for col in ["gene", "product", "description"]:
        if col in df.columns:
            df[col] = df[col].fillna("").astype(str)
    df["len"] = df["length"].astype(int)

    def is_n(row):
        s = f"{row['gene']} {row['product']} {row['description']}".lower()
        return any(k in s for k in ["nucleoprotein", "nucleocapsid", " protein n"])

    def is_l(row):
        s = f"{row['gene']} {row['product']} {row['description']}".lower()
        return any(k in s for k in ["rna-directed rna polymerase", " protein l ", "rd rp", "rdrp", "replicase", "transcriptase"]) or row['gene'].lower() == 'l'

    def is_gpc(row):
        s = f"{row['gene']} {row['product']} {row['description']}".lower()
        return any(k in s for k in ["glycoprotein precursor", "envelope polyprotein", "envelopment polyprotein", "m polyprotein", "glycoprotein (g1/g2)"]) or row['gene'].lower() in ['m','gpc']

    def is_gn(row):
        s = f"{row['gene']} {row['product']} {row['description']}".lower()
        return any(k in s for k in ["glycoprotein n", "gn", "g1"]) and "precursor" not in s

    def is_gc(row):
        s = f"{row['gene']} {row['product']} {row['description']}".lower()
        return any(k in s for k in ["glycoprotein c", "gc", "g2"]) and "precursor" not in s

    df["is_N"] = df.apply(is_n, axis=1)
    df["is_L"] = df.apply(is_l, axis=1)
    df["is_GPC"] = df.apply(is_gpc, axis=1)
    df["is_Gn_only"] = df.apply(is_gn, axis=1)
    df["is_Gc_only"] = df.apply(is_gc, axis=1)
    return df


def filter_lengths(df: pd.DataFrame, key: str) -> pd.DataFrame:
    if key == 'N':
        return df[(df['len'] >= 350) & (df['len'] <= 500)]
    if key == 'L':
        return df[(df['len'] >= 1700) & (df['len'] <= 2500)]
    if key == 'GPC':
        return df[(df['len'] >= 900) & (df['len'] <= 1400)]
    return df


def make_category_fasta(category: str, ncbi_seqs: Dict[str, str], df: pd.DataFrame, out_fa: str):
    # select rows
    if category == 'N':
        sel = df[df['is_N']]
    elif category == 'L':
        sel = df[df['is_L']]
    elif category == 'GPC':
        sel = df[df['is_GPC']]
    else:
        raise ValueError(category)
    sel = filter_lengths(sel, category)

    # deduplicate by sequence
    seen = set()
    lines = []
    # add UniProt reference first
    ref_acc = REFS[category]['acc']
    ref_header, ref_seq = get_uniprot_seq(ref_acc)
    lines.append(f">{ref_acc}|reference\n{ref_seq}\n")
    seen.add(ref_seq)

    # add NCBI sequences
    count = 0
    for acc in sel['accession']:
        seq = ncbi_seqs.get(acc)
        if not seq:
            continue
        if set(seq) - AASET - set(['-','X','*','?']):
            # contains non-standard AA, keep but strip unknowns
            seq = ''.join([a for a in seq if a in AASET])
        if len(seq) < 50:
            continue
        if seq in seen:
            continue
        seen.add(seq)
        lines.append(f">{acc}\n{seq}\n")
        count += 1
        if count >= 5000:
            break
    with open(out_fa, 'w') as f:
        f.writelines(lines)


def run_mafft(in_fa: str, out_fa: str):
    cmd = ["mafft", "--auto", "--anysymbol", in_fa]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    with open(out_fa, 'w') as f:
        f.write(res.stdout)


def compute_entropy_and_deltaH(aln_fa: str, ref_id_prefix: str) -> pd.DataFrame:
    # load alignment into a list of sequences
    seqs = []
    ref_seq = None
    ref_name = None
    for rec in SeqIO.parse(aln_fa, "fasta"):
        s = str(rec.seq).upper()
        seqs.append(s)
        if rec.id.startswith(ref_id_prefix):
            ref_seq = s
            ref_name = rec.id
    if ref_seq is None:
        # try matching by description
        for rec in SeqIO.parse(aln_fa, "fasta"):
            if ref_id_prefix in rec.description:
                ref_seq = str(rec.seq).upper()
                ref_name = rec.id
                break
    assert ref_seq is not None, f"Reference {ref_id_prefix} not found in alignment"

    L = len(seqs[0])
    for s in seqs:
        assert len(s) == L

    # compute per-column entropy H and deltaH = log2(20) - H
    H = []
    dH = []
    aa_cols = []
    for i in range(L):
        col = [s[i] for s in seqs if s[i] != '-']
        counts = Counter([a for a in col if a in AASET])
        n = sum(counts.values())
        if n == 0:
            H.append(0.0)
            dH.append(0.0)
        else:
            probs = [c / n for c in counts.values()]
            h = -sum(p * math.log2(p) for p in probs if p > 0)
            H.append(h)
            dH.append(LOG2_20 - h)
        aa_cols.append(ref_seq[i])

    # map alignment columns to reference positions (skip gaps in reference)
    rows = []
    ref_pos = 0
    for i, aa in enumerate(aa_cols):
        if aa == '-':
            continue
        ref_pos += 1
        rows.append({
            'pos': ref_pos,
            'residue': aa,
            'H': H[i],
            'deltaH': dH[i],
        })
    return pd.DataFrame(rows)


def main():
    ncbi_seqs = read_fasta_dict(NCBI_FASTA)
    df_meta = classify_ncbi_records()

    results = {}
    for cat in ['N', 'GPC', 'L']:
        out_fa = os.path.join(CONSDIR, f"{cat}.fasta")
        aln_fa = os.path.join(CONSDIR, f"{cat}.aln.fasta")
        make_category_fasta(cat, ncbi_seqs, df_meta, out_fa)
        run_mafft(out_fa, aln_fa)
        df_scores = compute_entropy_and_deltaH(aln_fa, REFS[cat]['acc'])
        out_csv = os.path.join(CONSDIR, f"{cat}_deltaH.csv")
        df_scores.to_csv(out_csv, index=False)
        results[cat] = out_csv

    # Derive Gn and Gc from GPC mapping using reference ranges
    gpc_csv = results['GPC']
    df_gpc = pd.read_csv(gpc_csv)
    gn_start, gn_end = REFS['GPC']['gn_range']
    gc_start, gc_end = REFS['GPC']['gc_range']
    df_gn = df_gpc[(df_gpc['pos'] >= gn_start) & (df_gpc['pos'] <= gn_end)].copy()
    df_gc = df_gpc[(df_gpc['pos'] >= gc_start) & (df_gpc['pos'] <= gc_end)].copy()
    df_gn['pos'] = df_gn['pos'] - gn_start + 1
    df_gc['pos'] = df_gc['pos'] - gc_start + 1
    gn_csv = os.path.join(CONSDIR, f"Gn_deltaH.csv")
    gc_csv = os.path.join(CONSDIR, f"Gc_deltaH.csv")
    df_gn.to_csv(gn_csv, index=False)
    df_gc.to_csv(gc_csv, index=False)

    manifest = {
        'N': results['N'],
        'GPC': results['GPC'],
        'Gn': gn_csv,
        'Gc': gc_csv,
        'L': results['L'],
    }
    with open(os.path.join(CONSDIR, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2)
    print(json.dumps(manifest, indent=2))

if __name__ == '__main__':
    main()
