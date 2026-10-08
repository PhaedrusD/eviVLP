# Case Study 2: Qβ with the Staphylococcus aureus Z Domain

Reproduces the hybrid VLP system described in Brown, Fiedler, and Finn, *Biochemistry* **48**(47), 11155–11157 (2009), in which the IgG-binding Z domain of *S. aureus* protein A is fused to the C-terminus of the Qβ bacteriophage coat protein via an eight-residue spacer. Hybrid capsids form when the fusion protein is co-expressed with the wild-type coat protein, displaying the Z domain on the capsid exterior.

## Inputs

| File | Description |
|---|---|
| `2SPZ.pdb` | NMR solution structure of the Z domain (RCSB PDB 2SPZ, model 1 extracted) |
| `2SPZ_P8.pdb` | Z domain with an N-terminal (G₄S)G₃ spacer; the structure used as the GUI input |

The Qβ coat protein template is not duplicated here; the GUI reads it from `../../genetic_modification/data/fix_a.pdb`, `fix_b.pdb`, and `fix_c.pdb` (the three quasi-equivalent T=3 conformations).

## Methodology

1. Obtain the Z domain starting structure from RCSB PDB entry **2SPZ**. 2SPZ contains 10 NMR conformers; extract model 1 (or the lowest-energy model) before further work.
2. Construct the fused sequence in FASTA format: Z domain sequence prepended (or appended, depending on construct orientation) by the chosen linker.
3. Predict the structure of the fused construct. Recommended tools:
   - **AlphaFold 3 server** ([alphafoldserver.com](https://alphafoldserver.com/)) — current state of the art, free tier
   - **ColabFold** ([github.com/sokrypton/ColabFold](https://github.com/sokrypton/ColabFold)) — runs AlphaFold 2 in Google Colab, no local install needed
   - **ESMFold** ([esmatlas.com](https://esmatlas.com/resources?action=fold)) — faster than AlphaFold for short sequences; lower accuracy on complex folds
   For short linkers attached to small, well-folded domains like the Z domain, any of these works. AlphaFold variants tend to predict flexible linkers as disordered, which is honest about the underlying physics.
4. Save the predicted structure as `2SPZ_P8.pdb` (or equivalent name) in this folder. Inspect the pLDDT scores around the linker — they will be low (~30–50), which is expected for a flexible peptide and not a problem for the downstream workflow.
5. Run the genetic modification GUI with `2SPZ_P8.pdb` as the protein insert and select Qβ as the capsid type. Set the number of fusion proteins to match the experimental condition being modeled (Brown et al. report ~6, ~20, or ~33 per particle depending on growth conditions).

## Scientific Record

Predictions in this case study were generated with AlphaFold2/ColabFold v1.6.1 on 6/11/26

## Limitations and assumptions

The original publication does not disclose the linker amino acid sequence beyond stating that it is an "octapeptide spacer." It is named CP8Z, perhaps it is **(Pro)₈** — a rigid, extended linker — based on the desire to project the Z domain away from the capsid surface for IgG binding. We preferred a flexible linker such as (G₄S)G₃ and it does not affect the downstream assembly workflow demonstrated here, only the predicted geometry of the unfused Z–linker–coat protein structure. Anyone reproducing this work with a definitive linker sequence (e.g., obtained from the corresponding author or the deposited pCCP8Z plasmid) can substitute that sequence at step 2 without other changes to the procedure. 

## Citation

Brown, S. D.; Fiedler, J. D.; Finn, M. G. Assembly of Hybrid Bacteriophage Qβ Virus-Like Particles. *Biochemistry* **2009**, *48* (47), 11155–11157. DOI: [10.1021/bi901306p](https://doi.org/10.1021/bi901306p). PMC: [PMC2799296](https://pmc.ncbi.nlm.nih.gov/articles/PMC2799296/).
