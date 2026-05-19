cat > ~/eviVLP-repo/README.md << 'EOF'
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
