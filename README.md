# Andes Orthohantavirus Endonuclease & Influenza Inhibitor Repurposing Pipeline

**Structural Conservation of the Andes Orthohantavirus L Endonuclease Active Site Enables Cross-Target Evaluation of Influenza Cap-Snatching Inhibitors**

This repository contains the computational workflow, processed datasets, structural analyses, molecular docking results, cheminformatics workflows, and supporting scripts associated with the study.

The project investigates the structural conservation of the endonuclease active site of Andes orthohantavirus (ANDV) and evaluates whether compounds previously associated with influenza A polymerase acidic (PA) endonuclease activity can be computationally assessed against the ANDV L-protein endonuclease.

---

## 1. Study Overview

Andes orthohantavirus (ANDV) is the causative agent of Hantavirus Pulmonary Syndrome (HPS), a severe respiratory disease associated with a reported case-fatality rate of approximately 30–50%. Currently, no approved small-molecule antiviral treatment specifically targets ANDV.

The N-terminal endonuclease domain of the viral L protein participates in the cap-snatching mechanism, cleaving host-derived RNA fragments to facilitate viral transcription. A related enzymatic mechanism is employed by the polymerase acidic (PA) subunit of influenza A virus.

This study investigates the structural similarities between these endonucleases and explores the potential cross-target compatibility of previously reported influenza PA-active compounds through computational analysis.

### Study Objectives

The computational workflow comprises four major components:

1. **Sequence Conservation Analysis:** Evaluation of Shannon entropy across ANDV structural and enzymatic proteins, including nucleoprotein (N), glycoprotein precursor (GPC), Gn, Gc, and RNA-dependent RNA polymerase (L).

2. **Comparative Structural Analysis:** Investigation of the binuclear manganese coordination environment and catalytic-site geometry across ANDV (PDB: `6Q99`), Hantaan virus (PDB: `5IZE`), and influenza A PA (PDB: `6FS6`).

3. **Cross-Target Molecular Docking:** Evaluation of a curated ChEMBL library of influenza PA-active compounds against both ANDV and influenza A endonuclease structures using an automated DockFlow/AutoDock Vina workflow.

4. **Cheminformatics and Analogue Generation:** Prioritization of parent scaffolds, physicochemical filtering, structural-alert assessment, and computational generation of analogues through fluorine substitution, methyl substitution, and ring-nitrogen replacement.

---

## 2. Key Results

### 2.1. ANDV Protein Sequence Conservation

Shannon entropy analysis was performed across multiple sequence alignments to characterize sequence variability among ANDV proteins.

| Protein | Alignment Depth | Positions | Mean ΔH | Median ΔH | Positions ≥ 4.0 (%) | Positions ≥ 3.5 (%) |
|:---|---:|---:|---:|---:|---:|---:|
| **N** | 27 | 428 | 4.24 | 4.32 | 91.36 | 96.03 |
| **GPC** | 61 | 1,138 | 4.25 | 4.32 | 91.65 | 98.24 |
| **Gn** | 61 | 633 | 4.23 | 4.32 | 89.10 | 97.79 |
| **Gc** | 61 | 487 | 4.28 | 4.32 | 95.28 | 98.97 |
| **L** | 6 | 2,153 | 4.20 | 4.32 | 86.44 | 94.43 |

**Note:** The L-protein alignment contained only six sequences, limiting the maximum observable Shannon entropy to log₂(6) = 2.585 bits. The entropy values and thresholds should therefore be interpreted in the context of the respective alignment depths and the exact normalization used in the analysis.

### 2.2. Comparative Structural Analysis and Metal Coordination

The binuclear manganese coordination environments of ANDV, influenza A PA, and Hantaan virus endonucleases were compared using structural superposition and geometric measurements.

| Structural Feature | ANDV (6Q99) | Influenza A (6FS6) | Hantaan Virus (5IZE) |
|:---|---:|---:|---:|
| Resolution (Å) | 2.95 | 2.29 | 1.70 |
| Mn–Mn Distance (Å) | 3.79 | 3.85 | 3.73 |
| Anchor Fit (10 atoms, Å) | — | 1.17 | 1.40 |
| Catalytic Core Fit (14 atoms, Å) | — | 1.47 | 1.22 |
| Matched Mn RMSD (Å) | — | 0.92 | 0.65 |

The structural comparisons indicate conservation of the binuclear metal arrangement and the overall geometry of the catalytic region across the examined endonucleases.

**Structural references:**
- ANDV endonuclease: PDB `6Q99`
- Influenza A PA endonuclease: PDB `6FS6`
- Hantaan virus endonuclease: PDB `5IZE`

### 2.3. Comparative Molecular Docking

A curated library of 246 influenza PA-active compounds was evaluated against both ANDV and influenza A endonuclease structures using AutoDock Vina.

| Docking Metric | ANDV (6Q99) | Influenza A (6FS6) |
|:---|---:|---:|
| Median Docking Score (kcal/mol) | −6.49 | −5.62 |
| Interquartile Range (kcal/mol) | −7.16 to −6.03 | −6.05 to −5.11 |
| Mean ± SD (kcal/mol) | −6.54 ± 0.85 | −5.59 ± 0.74 |
| Score Range (kcal/mol) | −8.26 to −3.93 | −7.60 to −3.56 |
| Spearman Correlation with Heavy-Atom Count (ρ) | −0.499 | −0.689* |

*The influenza A correlation value is formatted as −0.689 based on the supplied value of −689; verification against the original statistical output is required.*

The ANDV docking score distribution showed more negative median and mean values than the influenza A PA results. These differences represent computational docking-score comparisons and should not be interpreted as direct evidence of stronger experimental binding affinity.

### 2.4. Prioritized Parent Compounds

The following ten ChEMBL compounds represent a subset of the prioritized parent scaffolds identified through comparative docking and physicochemical profiling.

| Rank | ChEMBL ID | ANDV Score | PA Score | Δ Score | pChEMBL | MW (Da) | cLogP | TPSA (Å²) |
|---:|:---|---:|---:|---:|---:|---:|---:|---:|
| 1 | CHEMBL2040585 | −8.06 | −6.31 | −1.75 | 4.61 | 363 | 1.7 | 121 |
| 2 | CHEMBL2040582 | −8.01 | −6.65 | −1.36 | 5.79 | 351 | 1.7 | 121 |
| 3 | CHEMBL2040580 | −7.91 | −5.99 | −1.92 | 5.64 | 337 | 1.4 | 121 |
| 4 | CHEMBL2040581 | −7.91 | −6.00 | −1.91 | 6.54 | 337 | 1.4 | 121 |
| 5 | CHEMBL2040576 | −7.83 | −5.97 | −1.86 | 5.92 | 341 | 1.0 | 121 |
| 6 | CHEMBL2040578 | −7.83 | −6.15 | −1.67 | 6.14 | 337 | 1.4 | 121 |
| 7 | CHEMBL4438546 | −7.79 | −6.37 | −1.43 | 7.90 | 323 | 3.2 | 100 |
| 8 | CHEMBL2040574 | −7.79 | −5.95 | −1.84 | 5.78 | 323 | 0.9 | 121 |
| 9 | CHEMBL2040583 | −7.76 | −5.94 | −1.82 | 5.86 | 343 | 1.6 | 121 |
| 10 | CHEMBL2040575 | −7.71 | −6.11 | −1.60 | 5.15 | 337 | 1.2 | 121 |

All docking scores are reported in kcal/mol.

These compounds are computationally prioritized structures from the screened library, not experimentally validated ANDV inhibitors.

### 2.5. Structure-Based Analogue Generation

Three analogue-generation strategies were applied to the prioritized parent scaffolds.

| Modification Strategy | Enumerated | Novel | Passed Filters | Pass Rate (%) |
|:---|---:|---:|---:|---:|
| C–H to C–F substitution (Fluorine Scan) | 126 | 125 | 119 | 95.2 |
| C–H to C–CH₃ substitution (Methyl Scan) | 126 | 124 | 123 | 99.2 |
| Ring C–H to N substitution (Ring-Nitrogen Walk) | 126 | 126 | 70 | 55.6 |
| **Total** | **378** | **375** | **312** | **83.2** |

The generated structures were evaluated using computational novelty checks and physicochemical filtering criteria. The resulting analogue set provides a basis for further structure-based evaluation.

---

## 3. Computational Workflow

The project integrates sequence analysis, structural bioinformatics, molecular docking, and cheminformatics into a unified computational framework.

### 3.1. Sequence Conservation Analysis

- Collection of ANDV protein sequences from public sequence databases.
- Protein-specific sequence grouping and domain mapping.
- Multiple sequence alignment using MAFFT.
- Column-wise Shannon entropy calculations.
- Identification and visualization of conserved and variable positions.

Proteins investigated:
- Nucleoprotein (N)
- Glycoprotein precursor (GPC)
- Glycoprotein Gn
- Glycoprotein Gc
- RNA-dependent RNA polymerase (L)

### 3.2. Structural Superposition

Comparative structural analysis was performed using the following structures:

| Virus | Protein | PDB ID |
|:---|:---|:---|
| Andes orthohantavirus | L-protein endonuclease | `6Q99` |
| Influenza A virus | PA endonuclease | `6FS6` |
| Hantaan virus | L-protein endonuclease | `5IZE` |

The analysis included:

- Structural alignment using Kabsch-based superposition.
- Comparison of binuclear manganese coordination.
- Measurement of metal–metal distances.
- Evaluation of conserved metal-coordination anchors.
- Local catalytic-core RMSD calculations.
- Matched metal-position RMSD calculations.

### 3.3. Binding-Site and Pocket Analysis

Binding-site identification and characterization included:

- DeepSite-based pocket prediction.
- Pocket detection and characterization.
- Binding-site coordinate analysis.
- Comparison of predicted pockets with the catalytic region.
- Preparation of receptor structures and docking grids.

### 3.4. Molecular Docking

A curated ChEMBL compound library containing 246 paired docking results was evaluated against ANDV and influenza A PA endonuclease structures.

The automated docking workflow incorporated:

- Receptor preparation and structural preprocessing.
- PDBFixer-based structure preparation.
- Open Babel-based ligand and receptor format conversion.
- Docking-grid configuration.
- AutoDock Vina-based molecular docking.
- Extraction and comparison of docking scores.
- Statistical analysis of paired docking results.

The DockFlow workflow was used to organize and automate the docking process.

### 3.5. Cheminformatics and Physicochemical Filtering

Prioritized compounds were evaluated using molecular descriptors and structural-alert analysis.

The analysis included:

- Molecular weight (MW).
- Calculated lipophilicity (cLogP).
- Topological polar surface area (TPSA).
- Lipinski rule-of-five assessment.
- Structural-alert evaluation.
- Assessment of metal-binding and chelating functional groups.

Particular attention was given to the **chelator paradox**: metal-binding functional groups may contribute to interactions with the catalytic metal center but can also introduce nonspecific metal-binding, selectivity, and developability concerns.

### 3.6. Analogue Enumeration

Three structure-based modification strategies were implemented:

1. **Fluorine Scan:** Replacement of eligible carbon-bound hydrogen atoms with fluorine.
2. **Methyl Scan:** Replacement of eligible carbon-bound hydrogen atoms with methyl groups.
3. **Ring-Nitrogen Walk:** Replacement of eligible aromatic or ring carbon positions with nitrogen.

Generated analogues were subjected to novelty checks and physicochemical filtering.

---

## 4. Repository Structure

```text
.
├── analysis/                  # Processed results, joined datasets and statistics
├── data/                      # Protein sequences, alignments, PDB structures and entropy data
├── deepsite_6Q99/              # DeepSite pocket predictions and grid configurations
├── docking/                   # Docking configurations, receptors, grids and score outputs
├── docs/                      # Reports, domain mapping and pipeline documentation
├── ligands/                   # Ligand mappings and structural identifiers
├── pdb_cache/                 # Downloaded PDB/CIF reference structures
├── pharmacophore/              # Pharmacophore extraction and comparison
├── phylogenetics/              # Alignments, phylogenetic trees and evolutionary analyses
├── plots/                     # Generated figures and statistical visualizations
├── pockets/                   # Pocket detection and binding-site analyses
├── pocket_similarity/         # Pocket similarity rankings and filtering results
├── potential drugs/           # Prioritized compounds, analogues and filtering results
├── receptors/                 # Prepared receptor structures and coordinate information
├── reports/                   # Summary reports and statistical outputs
├── screenshots/               # Structural visualizations and pocket snapshots
├── scripts/                   # Python scripts and shell automation
├── tables/                    # Result tables, ADMET classifications and compound shortlists
└── materials_and_methods.txt  # Detailed computational protocol notes
```

*Note: The directory structure reflects the organization of the repository. Individual files and subdirectories may vary between versions.*

---

## 5. Software and Dependencies

The computational analyses were implemented using the following software and libraries:

| Software / Library | Application |
|:---|:---|
| Python 3.12 | Main programming environment |
| RDKit | Molecular processing and cheminformatics |
| Biopython | Protein sequence and structural bioinformatics |
| NumPy | Numerical computation |
| pandas | Dataset processing and analysis |
| SciPy | Statistical analysis |
| Matplotlib | Data visualization |
| MAFFT | Multiple sequence alignment |
| AutoDock Vina | Molecular docking |
| Open Babel | Molecular format conversion |
| PDBFixer | Protein structure preparation |
| DeepSite | Binding-pocket prediction |
| DockFlow | Docking workflow automation |

### Installation

Clone the repository:

```bash
git clone <REPOSITORY_URL>
cd <REPOSITORY_DIRECTORY>
```

Create and activate a Python virtual environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
```

Install the principal Python dependencies:

```bash
pip install numpy pandas scipy matplotlib biopython rdkit
```

Additional external programs, including MAFFT, AutoDock Vina, Open Babel, and the required DockFlow utilities, must be installed and configured separately.

**Note:** The installation commands provide the principal Python dependencies, not a guaranteed complete reproduction environment. Refer to the individual workflow scripts and `materials_and_methods.txt` for additional requirements and configuration details.

---

## 6. Data Sources and Structural References

The study relies on publicly available biological and chemical resources.

### Protein Structures

| PDB ID | Structure |
|:---|:---|
| `6Q99` | Andes orthohantavirus endonuclease |
| `6FS6` | Influenza A polymerase acidic subunit endonuclease |
| `5IZE` | Hantaan virus endonuclease |

Structural data are available through the [RCSB Protein Data Bank](https://www.rcsb.org/).

### Chemical Data

Compound information and activity annotations were obtained from [ChEMBL](https://www.ebi.ac.uk/chembl/).

---

## 7. Interpretation and Limitations

This repository documents a computational investigation of structural conservation and potential cross-target compatibility between ANDV and influenza A endonucleases.

The results should be interpreted within the following limitations:

- Sequence conservation does not independently establish functional equivalence.
- Structural similarity does not necessarily imply equivalent ligand-binding behavior.
- AutoDock Vina scores are approximate scoring-function outputs and are not experimentally measured binding affinities.
- Docking-score differences between targets do not establish antiviral activity or target selectivity.
- Metal-chelating interactions require particular caution because of potential nonspecific binding and limitations in conventional docking representations.
- Physicochemical filtering and analogue enumeration do not establish pharmacokinetic safety, toxicity, or biological efficacy.
- The sequence depth of the L-protein alignment was limited compared with the other analyzed proteins.

The prioritized compounds and generated analogues should therefore be considered **computational hypotheses for further investigation**, rather than validated antiviral candidates.

Further work would require experimental binding studies, biochemical endonuclease inhibition assays, antiviral activity testing, and appropriate toxicity evaluation.

---

## 8. Authors and Affiliations

**Aaditya Raj Sharma**¹ — Corresponding Author  
Email: aadityarajsharmamars@gmail.com

**Jayesh Mehra**¹

**Sana Khan**¹

¹ Invariant Genomics Private Limited, India.

---

## 9. Citation

If you use this repository or its computational workflow, please cite the associated study:

> Sharma, A. R., Mehra, J., & Khan, S. Structural conservation of the Andes orthohantavirus L endonuclease active site enables cross-target evaluation of influenza cap-snatching inhibitors.

*Publication details to be added upon publication.*

---

## 10. License and Availability

The protein structures used in this study are publicly available through the RCSB Protein Data Bank, and chemical compound information is accessible through ChEMBL.

The repository provides computational scripts, processed datasets, structural analyses, docking outputs, and supporting documentation for academic and research purposes.

A formal software license should be added to the repository to specify the terms under which the source code and associated materials may be reused, modified, and redistributed.

---

**Developed by Invariant Genomics Private Limited**  
*Computational Biology | Structural Bioinformatics | Cheminformatics | Antiviral Research*
