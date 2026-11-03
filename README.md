# eviVLP

A toolkit for visualizing and preparing chemically and genetically modified virus-like particles (VLPs) for molecular dynamics simulation.

## Overview

Virus-like particles are self-assembled nanoparticles derived from viral capsid proteins, with applications in imaging, drug delivery, vaccines, and gene therapy. Experimental work with VLPs often involves modifying coat proteins to display additional proteins, peptides, or small molecules — modifications typically targeted to specific solvent-exposed residues.

Accurately modeling such modified VLPs at atomic resolution is non-trivial. It requires reproducing covalent bonding patterns, respecting steric constraints from the full capsid architecture, and producing topology and parameter files compatible with the simulation force field. eviVLP automates these workflows, including the often-difficult task of reconciling parameters across multiple force fields (e.g., CHARMM36 protein types alongside CGenFF small-molecule types).

## Workflows

The repository contains two project workflows, each built around a wxWidgets GUI:

- **`chemical_modification/`** — attach small drug-like molecules or imaging agents to specific coat-protein residues (Lys, Cys, Tyr, N-terminus). Automates novel amino acid construction, topology and parameter file generation, and cross–force field parameter handling, with hooks to CGenFF, VMD/psfgen, and NAMD.
- **`genetic_modification/`** — fuse additional protein sequences onto coat proteins and assemble full VLP capsids (Qβ, TMV) with the modified subunit replicated across all symmetry-related copies.

## Repository structure

Each project contains:

- `src/` — C++ source (wxWidgets GUI and supporting code)
- `scripts/` — Python utilities invoked by the GUI
- `tcl/` — VMD / NAMD Tcl scripts
- `data/` — force-field files, structure templates and symmetry matrices (see "What is included" below)

## Quick start (chemical modification GUI)

```
git clone https://github.com/PhaedrusD/eviVLP.git
cd eviVLP
./setup.sh                 # run folder defaults to ~/eviVLP_run; pass a path to choose
cd ~/eviVLP_run && ./GUI
```

`setup.sh` checks the build prerequisites (installing them with Homebrew on macOS when it can), lays the repository's `src/`, `scripts/`, `tcl/` and `data/` out as the single flat run folder the GUI works in, builds the GUI and the internal-coordinate generators there, and reports where it found VMD and NAMD. Re-running it after a `git pull` refreshes the run folder without overwriting files you have edited there.

VMD and NAMD are located at startup — from `EVIVLP_VMD` / `EVIVLP_NAMD` if set, otherwise from the PATH, otherwise from the usual install locations (`/Applications/VMD*.app`, `~/Downloads/NAMD*`). The status log's first lines show which were chosen.

To update later: `git pull` in the clone, then `./setup.sh` again.

## Case studies

Two worked examples live in `case_studies/`, one per workflow. Each folder has its own `README.md` with the background, the step-by-step run, and notes on the choices made, plus an `inputs/` folder holding the files you load into the GUI.

```
cd case_studies/case_study1_TMV_TEMPO        # chemical modification: TEMPO spin label on TMV lysine
cat README.md                                # or open it in any editor / on GitHub
ls inputs/                                   # TEMPO.pdb is the drug; nitroxide_overrides.prm the charge-override file

cd ../case_study2_QB_zdomain                 # genetic modification: Z domain fused to the Qβ coat protein
cat README.md
ls inputs/                                   # 2SPZ_P8.pdb is the GUI input
```

If you are trying eviVLP for the first time, start with case study 1: it exercises the whole chemical-modification pipeline, and its `inputs/` folder is the test input referred to above. When the GUI asks for a file, browse to the `inputs/` folder of the case study; nothing needs to be copied into the run folder by hand.

## Requirements

Software you install:

- A C++17 compiler and wxWidgets 3.x (`brew install wxwidgets` on macOS)
- Python 3 with NumPy and Biopython (`pip install numpy biopython`); the genetic-modification GUI additionally needs MDAnalysis
- Open Babel (`brew install open-babel`) — converts the drug structure to MOL2
- VMD 1.9.4 or later (psfgen and viewing)
- NAMD 3.x (minimisation)
- An account on the [CGenFF web service](https://cgenff.silcsbio.com/) — the GUI writes the MOL2 files, you upload them there and import the returned stream files

`setup.sh` checks and, on macOS, installs the first three.

## What is included

You do not need to obtain any force-field files. Everything the pipelines read is in the repository and is laid out into the run folder by `setup.sh`:

| Module | File in `data/` | What it is |
|---|---|---|
| chemical | `top_all36_prot.rtf`, `par_all36m_prot.prm` | CHARMM36m protein topology and parameters (MacKerell lab release) |
| chemical | `par_water_ions.inp` | CHARMM TIP3P water and ion parameters |
| chemical | `par_all36_cgenff.prm` | CGenFF 4.6 parameters, with the NBFIX pair corrections against nucleic-acid types removed (`filter_cgenff_nbfix.py`) so the file loads alongside the protein set in NAMD |
| chemical | `top_all36_cgenff_CBD_template.rtf` | CGenFF 4.6 topology plus the eviVLP residues: `LYE`, `CYE`, `TYE` (modified Lys/Cys/Tyr scaffolds) and `BBQ`, `BBT` (the Qβ Ala1 and TMV Ser1 N-terminus models). **This is the modifiable topology.** Run Renaming regenerates the working `top_all36_cgenff_CBD.rtf` from it on every run, and Create Novel Amino Acid appends the run's `RESI NAA` / `RESI NPC` to that working copy — so edit the template, never the working copy |
| chemical | `par_GUI_template.prm` | seed for `par_GUI.prm`, the per-run parameter file that collects the CGenFF-derived and cross-force-field terms the released sets lack; regenerated by Run Renaming |
| chemical | `qb_cp{1,2,3}.pdb`, `qb_cp*_{cys,lys}.pdb`, `tmv_Y139.pdb`, `tmv_K53_K68.pdb` | coat-protein chain templates (Qβ A/B/C conformers; TMV with the Tyr139 or native Lys53/Lys68 sites marked `NAA`) |
| chemical | `matrix.txt`, `matrix_TMV.txt` | BIOMT symmetry operators (60 for the Qβ capsid, 49 for the TMV rod) |
| chemical | `patch_ZAP_*_template.txt`, `BBS_QB.pdb`, `conjugation_sites.dat` | patch templates per site type, the Qβ linker model, and the site registry |
| genetic | `top_all27_prot_na_CBD.inp` | combined CHARMM22/27 topology with the fusion residues the genetic GUI uses |
| genetic | `TMV_template.pdb`, `fix_{a,b,c}.pdb`, `matrix.txt`, `matrix_TMV.txt` | coat-protein templates and symmetry operators |

Files that are produced during a run — `one.mol2`, `naa.str`, `ntrm.str`, `nad.str`, `ntrm_clean.str`, `par_GUI.prm`, `top_all36_cgenff_CBD.rtf`, the filled `patch_maker_*.py` / `orient_place_drug_*.py` scripts, and every PSF/PDB the pipeline builds — are not tracked. The one input you supply is the drug or label structure as a PDB with all hydrogens; the case studies each provide one.

The CHARMM and CGenFF files are redistributed unmodified apart from the additions noted above; they are the work of the MacKerell laboratory (charmm.org, cgenff.com) and remain subject to their terms.

## CGenFF parameterization

Small molecules added to a VLP are parameterized using the CHARMM General Force Field (CGenFF), the CHARMM family's parameter set for non-polymer organic molecules. CGenFF assigns atom types, partial charges, and bonded/non-bonded parameters from its primary library — derived from QM and condensed-phase data and continuously expanded as new fragments are parameterized — and from analogy-based assignments (with penalty scores) where direct parameters are not yet available. The CGenFF program reads a MOL2-format structure and produces a CHARMM-compatible stream (STR) file containing the residue topology together with any parameters not already in the current CGenFF release. 

eviVLP automates the structural preparation and bookkeeping around this step. After the user supplies the modifying small molecule, the program converts it to MOL2 format with proper atom typing, then proceeds along one of two paths depending on the host operating system:

- **Linux:** the CGenFF command-line application is invoked directly, generating the STR file in place with no user intervention. The resulting topology and parameters are then integrated into the novel-amino-acid construction and merged into a custom parameter file alongside CHARMM36 protein terms.
- **macOS / Windows:** the CGenFF command-line application is not currently available on these platforms. The GUI instead writes the MOL2 file to disk and pauses, prompting the user to upload it to the SilcsBio CGenFF web service ([cgenff.silcsbio.com](https://cgenff.silcsbio.com/)), download the returned STR file, and place it in the working directory. The GUI then resumes the same downstream integration steps as the Linux path.

Both paths produce equivalent topology and parameter outputs. The Linux path is faster and fully automated; the web path requires a brief manual step but does not require local installation of CGenFF.

## Status

In active development as part of a PhD dissertation. A manuscript describing the software and its application to VLP–small molecule conjugates is in preparation.

## Author

**Cary Darwin**  
PhD Candidate, Nielsen Lab  
Department of Chemistry, University of Texas at Dallas

## Citation

Paper forthcoming
