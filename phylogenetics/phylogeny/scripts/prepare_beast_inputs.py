from pathlib import Path
import re
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

base = Path('phylogeny_expanded_v2')
ali_fa = base / 'rdRp_aligned_expanded.fasta'
meta_tsv = base / 'metadata.tsv'
outdir = Path('beast2_inputs')
outdir.mkdir(exist_ok=True, parents=True)

# Read alignment and set molecule type to protein for NEXUS export
aln = AlignIO.read(str(ali_fa), 'fasta')
prot_records = []
for rec in aln:
    r = SeqRecord(Seq(str(rec.seq)), id=rec.id, description="")
    r.annotations["molecule_type"] = "protein"
    prot_records.append(r)
prot_aln = MultipleSeqAlignment(prot_records)

AlignIO.write(prot_aln, str(outdir / 'rdRp_aligned.nexus'), 'nexus')

# Parse metadata to collect tip dates
label_to_year = {}
if meta_tsv.exists():
    lines = meta_tsv.read_text().splitlines()
    for line in lines[1:]:
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        label, header, year = parts[0], parts[1], parts[2]
        year = year.strip()
        if year:
            try:
                y = int(year)
                label_to_year[label] = y
            except ValueError:
                pass

# Write tip dates TSV
with (outdir / 'tipdates.tsv').open('w') as fh:
    fh.write('taxon\tdecimal_date\n')
    for rec in prot_aln:
        lab = rec.id
        if lab in label_to_year:
            fh.write(f"{lab}\t{float(label_to_year[lab])}\n")

# Write suggested calibrations (arthropod host-based, very broad priors)
calib_md = outdir / 'calibrations.md'
calib_md.write_text(
"""
Arthropod host-based calibration suggestions (broad, soft priors)

IMPORTANT: These priors are heuristic and reflect minimum ages based on arthropod host fossil records and major lineage splits. They are intentionally broad to avoid overconstraining rates. Review and adjust before analysis.

Suggested clades (use Monophyly + MRCA prior in BEAUti > Priors):

1) Tick-borne orthomyxoviruses (Thogotovirus + Quaranfil virus)
   - Taxa: Thogoto_virus_PB1, Quaranfil_virus_PB1
   - Rationale: Both primarily tick-associated. Oldest tick fossils ~90–100 Ma; Argasidae–Ixodidae split older (~200–300+ Ma).
   - Prior: LogNormal(meanInRealSpace=true, M=60, S=1.0, offset=90) Ma  [median ~60 Ma above a 90 Ma offset]

2) Tick-borne bunyavirus lineage (Nairoviridae) vs insect-borne bunyaviruses
   - Taxa (tick): Crimean_Congo_hemorrhagic_fever_virus_L
   - Taxa (insect/others): Bunyamwera_virus_L, Rift_Valley_fever_virus_L, Tomato_spotted_wilt_virus_L, Puumala_virus_L, Hantaan_virus_L, Andes_orthohantavirus_L
   - Apply MRCA prior to the node uniting all bunyaviruses (or specifically the split leading to Nairoviridae), with a soft minimum around the chelicerate–insect divergence.
   - Prior: LogNormal(meanInRealSpace=true, M=150, S=1.0, offset=250) Ma

3) Plant/thrips-associated tospoviruses vs mosquito-borne peribunya/phenu lineage
   - Taxa (tospovirus): Tomato_spotted_wilt_virus_L
   - Taxa (mosquito-borne): Bunyamwera_virus_L, Rift_Valley_fever_virus_L
   - Prior: LogNormal(meanInRealSpace=true, M=80, S=1.0, offset=100) Ma

Global model recommendations
- Substitution model: WAG + Gamma(4) + empirical AA frequencies
- Clock model: Uncorrelated Log-Normal (UCLN)
- Tree prior: Yule or Birth–Death (use Yule if uncertain)
- Tip dates: use tipdates.tsv for sequences with known sampling years; leave others undated.

Workflow in BEAUti (BEAST 2):
1. Import rdRp_aligned.nexus (Data type: Amino Acids).
2. Set clock to Relaxed Clock Log Normal.
3. Set Site Model to WAG + Gamma (4 categories) + Estimated frequencies.
4. Define Taxon Sets for the clades above (Monophyletic = true).
5. In Priors, add MRCA priors for each Taxon Set and configure LogNormal priors with the parameters above (units in millions of years). Use very wide S to keep constraints soft.
6. In Tip Dates, load tipdates.tsv (leave others blank).
7. Set chain length based on resources (e.g., 5e7 states), store every 5e3, burn-in 10–20%.
8. Run BEAST 2, then assess convergence in Tracer (ESS > 200) and summarize trees with TreeAnnotator.

Caveats
- Co-divergence with arthropod hosts is debated for RNA viruses; these priors may push rates unrealistically low. Inspect model fit and consider alternative calibrations (e.g., fossil-calibrated secondary priors or within-order calibrations) and partitioning by order.
- Consider running separate analyses for Orthomyxoviridae and Bunyavirales, then combining via secondary calibrations to mitigate time-dependent rate effects.
"""
)

(outdir / 'README.txt').write_text(
"""
BEAST 2 input package prepared.

Files:
- rdRp_aligned.nexus: Amino-acid alignment for PB1/L polymerases.
- tipdates.tsv: Sampling years for a subset of taxa (decimal years).
- calibrations.md: Suggested arthropod host-based node priors and workflow.

Open BEAUti (BEAST 2 GUI) and follow calibrations.md to set up the analysis.
"""
)

print('Wrote BEAST2 inputs to', outdir)
