#!/usr/bin/env python3
import os
import math
from typing import Dict, Tuple
import pandas as pd
from Bio.PDB import PDBParser, PDBIO
from Bio import pairwise2

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
DATADIR = os.path.join(BASE, "data")
CONSDIR = os.path.join(DATADIR, "conservation")
PDBDIR = os.path.join(DATADIR, "pdb")
OUTDIR = os.path.join(DATADIR, "pdb_colored")
os.makedirs(OUTDIR, exist_ok=True)

# Load per-protein deltaH tables
DELTAH = {
    'N': pd.read_csv(os.path.join(CONSDIR, 'N_deltaH.csv')).set_index('pos'),
    'GPC': pd.read_csv(os.path.join(CONSDIR, 'GPC_deltaH.csv')).set_index('pos'),
    'L': pd.read_csv(os.path.join(CONSDIR, 'L_deltaH.csv')).set_index('pos'),
}
MAX_DH = 4.321928094887363  # log2(20)

# Reference FASTA sequences (from our previous downloaded UniProt FASTA)
UNIPROT_FASTA = os.path.join(DATADIR, 'uniprot', 'andes_orthohantavirus_uniprot.fasta')
REFSEQ = {}
from Bio import SeqIO
for rec in SeqIO.parse(UNIPROT_FASTA, 'fasta'):
    if '|O36307|' in rec.description:
        REFSEQ['N'] = str(rec.seq).upper()
    if '|Q9E006|' in rec.description:
        REFSEQ['GPC'] = str(rec.seq).upper()
    if '|Q9E005|' in rec.description:
        REFSEQ['L'] = str(rec.seq).upper()

# Map PDB files to protein category
PDB_MAP = {
    '5E04_assembly1.pdb': 'N',
    '2K48_assembly1.pdb': 'N',
    '6Y5F_assembly1.pdb': 'GPC',
    '6Y5W_assembly1.pdb': 'GPC',
    '6Y6Q_assembly1.pdb': 'GPC',
    '6YRB_assembly1.pdb': 'GPC',
    '6YRQ_assembly1.pdb': 'GPC',
    '6ZJM_assembly1.pdb': 'GPC',
    '8DBZ_assembly1.pdb': 'GPC',
    '6Q99_assembly1.pdb': 'L',
}

AASET = set('ACDEFGHIKLMNPQRSWTYV')

def seq_from_chain(chain) -> str:
    seq = []
    three_to_one = {
        'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
        'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
    }
    for res in chain.get_residues():
        if res.get_id()[0] != ' ':
            continue
        name = res.get_resname()
        aa = three_to_one.get(name, 'X')
        seq.append(aa)
    return ''.join(seq)


def map_chain_positions_to_uniprot(chain_seq: str, ref_seq: str) -> Dict[int, int]:
    # global alignment with penalties allowing gaps
    aln = pairwise2.align.globalms(chain_seq, ref_seq, 2, -1, -5, -0.5, one_alignment_only=True)[0]
    a_chain, a_ref, score, start, end = aln
    mapping = {}
    cpos = 0
    rpos = 0
    for i in range(len(a_chain)):
        c = a_chain[i]
        r = a_ref[i]
        if c != '-':
            cpos += 1
        if r != '-':
            rpos += 1
        if c != '-' and r != '-':
            mapping[cpos] = rpos
    return mapping


def set_bfactor_by_variability(structure, chain_map: Dict[str, Dict[int, float]]):
    for model in structure:
        for chain in model:
            cid = chain.id
            var_by_pos = chain_map.get(cid, {})
            cpos = 0
            for res in chain:
                if res.get_id()[0] != ' ':
                    continue
                cpos += 1
                var = var_by_pos.get(cpos)
                for atom in res:
                    if var is None:
                        atom.set_bfactor(0.0)
                    else:
                        atom.set_bfactor(float(var))


def process_pdb(infile: str, protein: str) -> str:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('x', infile)

    ref_seq = REFSEQ[protein]
    dh = DELTAH[protein]

    chain_maps = {}
    # build mapping and compute variability per chain
    for model in structure:
        for chain in model:
            chain_seq = seq_from_chain(chain)
            if len(chain_seq) < 20:
                continue
            c2r = map_chain_positions_to_uniprot(chain_seq, ref_seq)
            varpos = {}
            for cpos, rpos in c2r.items():
                if rpos in dh.index:
                    deltaH = dh.loc[rpos, 'deltaH']
                    # variability = 1 - normalized conservation
                    var = max(0.0, 1.0 - (float(deltaH) / MAX_DH))
                    varpos[cpos] = var
            chain_maps[chain.id] = varpos
    set_bfactor_by_variability(structure, chain_maps)

    outfile = os.path.join(OUTDIR, os.path.basename(infile).replace('.pdb', '_var.pdb'))
    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile)
    return outfile


def main():
    outfiles = []
    for fname, prot in PDB_MAP.items():
        path = os.path.join(PDBDIR, fname)
        if not os.path.exists(path):
            continue
        out = process_pdb(path, prot)
        outfiles.append(out)
    print("\n".join(outfiles))

if __name__ == '__main__':
    main()
