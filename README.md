# eviVLP

A toolkit for visualizing and preparing chemically and genetically modified virus-like particles (VLPs) for molecular dynamics simulation.

## Overview

Virus-like particles are self-assembled nanoparticles derived from viral capsid proteins, with applications in imaging, drug delivery, vaccines, and gene therapy. Experimental work with VLPs often involves modifying coat proteins to display additional proteins, peptides, or small molecules — modifications typically targeted to specific solvent-exposed residues.

Accurately modeling such modified VLPs at atomic resolution is non-trivial. It requires reproducing covalent bonding patterns, respecting steric constraints from the full capsid architecture, and producing topology and parameter files compatible with the simulation force field. eviVLP automates these workflows, including the often-difficult task of reconciling parameters across multiple force fields (e.g., CHARMM36 protein types alongside CGenFF small-molecule types).

## Workflows

The repository contains two project workflows, each built around a wxWidgets GUI:

- **`chemical_modification/`** — attach small drug-like molecules or imaging agents to specific coat-protein residues (Lys, Cys, Tyr, N-terminus). Automates novel amino acid construction, topology and parameter file generation, and cross–force field parameter handling, with hooks to CGenFF, VMD/psfgen, and NAMD.
- **`genetic_modification/`** — fuse additional protein sequences onto coat proteins and assemble full VLP capsids (Qβ, TMV) with the modified subunit replicated across all symmetry-related copies.

## CGenFF parameterization

Small molecules added to a VLP are parameterized using the CHARMM General Force Field (CGenFF), which assigns atom types, partial charges, and bonded/non-bonded parameters by analogy to chemically similar fragments in its training set. The CGenFF program reads a MOL2-format structure and produces a CHARMM-compatible stream (STR) file containing the residue topology and any parameters not already present in the standard CGenFF release.

eviVLP automates the structural preparation and bookkeeping around this step. After the user supplies the modifying small molecule, the program converts it to MOL2 format with proper atom typing, then proceeds along one of two paths depending on the host operating system:

- **Linux:** the CGenFF command-line application is invoked directly, generating the STR file in place with no user intervention. The resulting topology and parameters are then integrated into the novel-amino-acid construction and merged into a custom parameter file alongside CHARMM36 protein terms.
- **macOS / Windows:** the CGenFF command-line application is not currently available on these platforms. The GUI instead writes the MOL2 file to disk and pauses, prompting the user to upload it to the SilcsBio CGenFF web service ([cgenff.silcsbio.com](https://cgenff.silcsbio.com/)), download the returned STR file, and place it in the working directory. The GUI then resumes the same downstream integration steps as the Linux path.

Both paths produce equivalent topology and parameter outputs. The Linux path is faster and fully automated; the web path requires a brief manual step but does not require local installation of CGenFF.

## Repository structure

Each project contains:

- `src/` — C++ source (wxWidgets GUI and supporting code)
- `scripts/` — Python utilities invoked by the GUI
- `tcl/` — VMD / NAMD Tcl scripts
- `data/` — force-field templates and topology files (chemical project only)

## Requirements

- wxWidgets 3.x
- Python 3 with NumPy
- VMD 1.9.4 or later
- NAMD 3.x
- CHARMM36 protein topology and parameters; CGenFF for small molecules
- CGenFF access — either the [SilcsBio web service](https://cgenff.silcsbio.com/) (macOS workflow) or the local command-line application (Linux workflow; not yet integrated in current build)

## Status

In active development as part of a PhD dissertation. A manuscript describing the software and its application to VLP–small molecule conjugates is in preparation.

## Author

**Cary Darwin**  
PhD Candidate, Nielsen Lab  
Department of Chemistry, University of Texas at Dallas

## Citation

Paper forthcoming
