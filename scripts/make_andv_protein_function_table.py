#!/usr/bin/env python3
import os
import re
import json
import textwrap
import requests
import pandas as pd

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
OUTDIR = os.path.join(BASE, "data", "summary")
os.makedirs(OUTDIR, exist_ok=True)

USER_AGENT = "PlayMoleculeAI/1.0 (Acellera; data-retrieval; contact: support@acellera.com)"

def get_json(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    r = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=60)
    r.raise_for_status()
    return r.json()


def extract_function(ujson):
    funs = []
    for c in ujson.get("comments", []):
        if c.get("commentType") == "FUNCTION":
            for tx in c.get("texts", []):
                v = tx.get("value")
                if v:
                    funs.append(v)
    txt = " ".join(funs)
    # collapse whitespace
    txt = re.sub(r"\s+", " ", txt).strip()
    # Keep first ~400 chars for readability
    if len(txt) > 600:
        txt = txt[:600].rsplit(". ", 1)[0] + "."
    return txt


def extract_chains(ujson):
    chains = []
    for feat in ujson.get("features", []):
        if str(feat.get("type")).lower() == "chain":
            desc = feat.get("description")
            loc = feat.get("location", {})
            begin = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            if begin and end:
                chains.append({
                    "description": desc,
                    "start": int(begin),
                    "end": int(end),
                })
    return chains


def primary_name(ujson):
    try:
        return ujson["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except Exception:
        # fallback to submitted name or first proteinName
        for k in ("submittedName", "alternativeNames"):
            if k in ujson.get("proteinDescription", {}):
                block = ujson["proteinDescription"][k]
                if isinstance(block, list) and block:
                    full = block[0].get("fullName", {}).get("value")
                    if full:
                        return full
        return None


def length_of_seq(ujson):
    try:
        return int(ujson["sequence"]["length"])
    except Exception:
        return None


def main():
    # Reviewed accessions for Andes orthohantavirus (from UniProt TSV):
    # O36307 (N), Q9E006 (GPC), Q9E005 (L)
    entries = [
        ("O36307", "S"),  # Nucleoprotein on S segment
        ("Q9E006", "M"),  # Glycoprotein precursor on M segment
        ("Q9E005", "L"),  # RdRp on L segment
    ]
    rows = []
    for acc, segment in entries:
        uj = get_json(acc)
        pname = primary_name(uj) or ""
        func = extract_function(uj)
        length = length_of_seq(uj)
        chains = extract_chains(uj)
        if acc == "Q9E006" and chains:
            # GPC: split into Gn/Gc chains (skip the overarching 'Envelopment polyprotein')
            for ch in chains:
                desc = (ch.get("description") or "").strip()
                if desc.lower().startswith("envelopment polyprotein"):
                    continue
                rows.append({
                    "protein": desc if desc else pname,
                    "segment": segment,
                    "uniprot": acc,
                    "positions": f"{ch['start']}-{ch['end']}",
                    "length": ch['end'] - ch['start'] + 1,
                    "function": func,
                })
        else:
            # Single-chain proteins (N, L)
            pos = None
            # Prefer chain feature if present
            if chains:
                ch0 = [c for c in chains if c.get("description") and c["description"].lower().startswith(pname.split()[0].lower())]
                ch = ch0[0] if ch0 else chains[0]
                pos = f"{ch['start']}-{ch['end']}"
                length = ch['end'] - ch['start'] + 1
            else:
                pos = f"1-{length}" if length else ""
            rows.append({
                "protein": pname,
                "segment": segment,
                "uniprot": acc,
                "positions": pos,
                "length": length,
                "function": func,
            })
    df = pd.DataFrame(rows, columns=["protein", "segment", "uniprot", "positions", "length", "function"])
    out_csv = os.path.join(OUTDIR, "andv_proteins_function_positions.csv")
    df.to_csv(out_csv, index=False)
    print(out_csv)

if __name__ == "__main__":
    main()
