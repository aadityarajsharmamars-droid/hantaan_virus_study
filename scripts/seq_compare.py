from Bio.PDB import PDBParser, PPBuilder
from Bio import Align
import sys

pdb1=sys.argv[1]
pdb2=sys.argv[2]

ppb=PPBuilder()
parser=PDBParser(QUIET=True)

s1=''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(parser.get_structure('A', pdb1))])
s2=''.join([str(pp.get_sequence()) for pp in ppb.build_peptides(parser.get_structure('B', pdb2))])

aligner=Align.PairwiseAligner()
aligner.mode='global'
aligner.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')
aligner.open_gap_score=-10
aligner.extend_gap_score=-0.5

alns=aligner.align(s1, s2)
if alns:
    a=alns[0]
    aln1=a.aligned[0]
    aln2=a.aligned[1]
    # compute identity over aligned positions
    i=j=0
    ident=0
    length=0
    for (s1b,s1e),(s2b,s2e) in zip(aln1, aln2):
        for k in range(min(s1e-s1b, s2e-s2b)):
            if s1[s1b+k]==s2[s2b+k]:
                ident+=1
            length+=1
    pid=100.0*ident/length if length>0 else 0.0
    print('len1=',len(s1),'len2=',len(s2))
    print('aligned_length=',length,'pid=%.2f%%'%pid,'score=',a.score)
else:
    print('No alignment produced')
