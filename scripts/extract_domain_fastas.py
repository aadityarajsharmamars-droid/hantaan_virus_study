from Bio.PDB import PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys, os

pairs=[('6FS6','receptors/6FS6/output.pdb'),('6Q99','receptors/6Q99/output.pdb')]
outdir='analysis/sequences'
os.makedirs(outdir, exist_ok=True)
parser=PDBParser(QUIET=True)
ppb=PPBuilder()

for tag, path in pairs:
    st=parser.get_structure(tag, path)
    # concatenate peptides per chain, keep first chain
    chains=list(st.get_chains())
    if not chains:
        continue
    # build peptides per chain and pick the longest chain
    best_seq=''; best_chain_id=None
    for ch in st.get_chains():
        seqs=[str(pp.get_sequence()) for pp in ppb.build_peptides(ch)]
        seq=''.join(seqs)
        if len(seq)>len(best_seq):
            best_seq=seq; best_chain_id=ch.id
    rec=SeqRecord(Seq(best_seq), id=f'{tag}_{best_chain_id}', description=f'{tag} chain {best_chain_id} (PDB-derived)')
    fn=os.path.join(outdir, f'{tag}_chain{best_chain_id}.fasta')
    SeqIO.write([rec], fn, 'fasta')
    print('Wrote', fn, 'len', len(best_seq))
