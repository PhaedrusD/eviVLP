# Case Study 1: kTMV with TEMPO (T158K, EDC coupling)

Reproduces the **kTMV–TEMPO** construct of Lumata *et al.*, in which the spin
label TEMPO is attached to a lysine introduced at the C-terminus of the TMV coat
protein by the single point mutation **T158K**. The nitroxide serves as an MRI
contrast agent that is quenched/reduced in the presence of superoxide.

The same publication also reports **nTMV–TEMPO**, in which TEMPO-azide is
attached to wild-type **Tyr139** by diazonium coupling followed by CuAAC. That
arm is not covered here. The two sites are 24.8 Å apart, and the paper
attributes kTMV's improved performance to the greater solvent exposure of the
lysine site.

## Inputs

In `inputs/`:

| File | Contents |
|---|---|
| `TEMPO.pdb` | 4-amino-TEMPO(H+) with a 6-aminohexanoic acid linker, 50 atoms, all hydrogens; linker amine protonated, so the conjugated residue is +1 |
| `fix_tempo_input.py` | how `TEMPO.pdb` was repaired; see "Repairing the input structure" |
| `nitroxide_overrides.prm` | hand-curated aminoxyl radical parameters; load after `par_GUI.prm` |
| `fix_atom_naming.py` | reconciles PDB atom names with the stream file |
| `tmv_K158.pdb` | the T158K coat protein; see "The mutated site" below |
| `TEMPO.mol2` | **stale** — 43-atom upload, predates the methyl fix |
| `TEMPO.str` | **stale** — CGenFF result for the 43-atom structure |
| `topology_entry.rtf` | **stale** — `RESI ONE` block from that same result |

The three files marked stale were generated from the 2,6-dimethyl structure.
They are kept only as a record of the earlier run. Re-running the CGenFF step
on the corrected `TEMPO.pdb` regenerates all three, and the atom names will
differ because the molecule now has six more atoms.

### The mutated site, and why it is not the GUI default

`inputs/tmv_K158.pdb` is `tmv_resall.pdb` with residue 158 renamed `NAA`;
threonine's `OG1` and `CG2` are removed, leaving `N, CA, C, CB, O` shared with
lysine. The lysine side chain (`CG, CD, CE, NZ`) and the whole TEMPO adduct are
then built by psfgen's `guesscoord` from the `RESI NAA` internal-coordinate
table, so the side-chain conformation is a model rather than a crystallographic
observation and must be minimised before use.

kTMV is a **T158K single-point mutant** -- wild-type TMV has no lysine at 158.
That makes it the wrong thing for the base GUI to ship as its lysine example,
for four separate reasons:

* it is one specific published construct, not a general capability;
* it requires a mutation, so the template cannot be derived from the wild-type
  structure by renaming alone;
* residue 158 is the C-terminus, so psfgen applies `CTER` to the `NAA` residue
  itself -- a path nothing else exercises;
* the entire side chain has to be guessed, which is the modelled-rotamer
  caveat above.

So the GUI's lysine default is now `tmv_K53_K68.pdb` -- the **native** Lys53 and
Lys68 of wild-type TMV (the doxorubicin sites of Finbloom et al. 2016). Those
keep their crystallographic `CG/CD/CE/NZ`, so `guesscoord` builds only the drug,
and residue 158 stays `THR` so `CTER` lands on an ordinary residue. None of the
four caveats above apply to it.

**To reproduce this case study**, copy `inputs/tmv_K158.pdb` beside the GUI
executable and change `GetTmvCoatProtein` in `chemical_modification/src/
FileProcessor.cpp` to return `"tmv_K158.pdb"` for `lys`. The comment there says
the same. The modelled-rotamer and `CTER`-on-`NAA` caveats apply only to this
arm, not to a default run.

### Relevant atoms of `RESI ONE`

From the **stale** `TEMPO.str`. The atom picks below will change when
CGenFF is re-run on the corrected structure; re-take them before use.

| Group | Atoms |
|---|---|
| Nitroxide (the reporter) | `N1` (NG301) – `O1` (OG312), on the piperidine ring `C1 C2 C3 N1 C5 C6` |
| Carboxylic acid (the handle) | `C13` (CG2O2) = `O3` (OG2D1), `C13` – `O2` (OG311) – `H25` |
| Link to the ring | `C12` – `C13` |

## Methodology

EDC activates the TEMPO carboxylic acid, which then acylates the lysine
ε-amine to give an amide. This is the acylation the `ZAP` patch models, so TEMPO
is supplied as the **free acid**, not as a pre-activated ester.

**Prerequisite (one time).** The TMV N-terminus model residue must be installed,
because a lysine job also modifies the N-terminus (below). Send
`chemical_modification/data/BBS_TMV.mol2` to CGenFF once, then
`python3 install_nterm_reference.py BBT.str` and click *Run Renaming*.

1. Select **TMV** as the VLP and **Lysine** as the residue. The GUI resolves the
   coat protein to `tmv_K158.pdb`.
2. Load `TEMPO.pdb` as the incoming molecule, and send the generated `one.mol2`
   to CGenFF. Import the returned zip; the GUI extracts `one.str` and verifies it
   against `one.mol2`.
3. Select the junction atoms in the viewer:

   | Slot | Atom | Role |
   |---|---|---|
   | 1st bonded (`*T1`) | `C13` | carbonyl carbon, bonds to `NZ` |
   | 2nd bonded (`*T2`) | `O3` | carbonyl oxygen |
   | 3rd bonded (`*T3`) | `C12` | α-carbon |
   | deleted | `O2`, `H25` | the hydroxyl leaving group |

   All three bonded slots must be filled; with fewer, the deleted-atom slots
   shift.
4. *Process File* orients TEMPO on the lysine (and, in parallel, on the
   N-terminus model).
5. *Process Patch Maker* writes `PRES ZAP` and `PRES ZPP` into the topology file
   and runs psfgen to give `patched.pdb/psf` and `ntrm.pdb/psf`.
6. Send **both** `naa.mol2` and `ntrm.mol2` to CGenFF — two separate jobs — and
   import both results. The GUI rejects a stream file that does not match the
   mol2 it was requested for.
7. *Create Novel Amino Acid* builds `RESI NAA` and `RESI NPC`, inserts the
   cross-force-field parameters, and validates both internal-coordinate tables.
8. *Create And View Capsid* builds the modified monomer, replicates it along the
   helical axis (49 copies), assembles the rod, then solvates and ionises it.

### Two consequences specific to this site

**The N-terminus is modified as well.** EDC coupling cannot distinguish the
lysine ε-amine from the free N-terminal α-amine, so the run also produces
`RESI NPC` at residue 1. TMV residue 1 is **serine**, so the N-terminus model is
the serine variant (`BBS_TMV.pdb`, `RESI BBT`, rename map `ntrm_ser`) rather
than the alanine model used for Qβ. The α-amine is free because the coat protein
is expressed recombinantly in *E. coli*, which lacks the eukaryotic
NatA/NatB/NatC N-terminal acetyltransferases; a plant-derived capsid may be
N-terminally acetylated and therefore unreactive.

**Residue 158 is the C-terminus.** psfgen applies its default `CTER` patch to the
`NAA` residue, which deletes `O`, adds `OT1`/`OT2`, and retypes `C` as `CC`.
`RESI NAA` provides the `CA`, `C` and `O` that `CTER` requires. Every other NAA
built by this toolkit has been an internal residue, so this is the first
terminal one.

## The radical oxygen

TEMPO's nitroxide `>N-O` is the reason the molecule is useful and the one group
CGenFF cannot type. No CGenFF release (v4.6 or v5.0 checked) has an aminoxyl
radical type. Left to itself, the CGenFF program reads the `N(+)-O(-)`
zwitterion that Open Babel writes, types the oxygen as `OG312`
("anionic alcohol oxygen"), and hands back a charge near `-0.66` on the oxygen
plus a spurious net `-1` on the whole residue. CGenFF v5.0 does now carry
`NG301 OG311`, but as an N-hydroxylamine: `1.4104 A` with an `N-O-H` angle
term. The aminoxyl bond is `~1.26 A` and has no hydrogen.

So the group is parameterised by hand. The values below are the ones actually
used in the 2023 production runs, recovered from `extra.prm` and
`RESI LYT` / `RESI YMT` in `top_all27_prot_na_CBD.inp`, and confirmed against
the production structure `T158K_new.psf` (segment `PA1`, residue 158). The
T158K and Y139 arms used identical values.

| | type | charge |
|---|---|---|
| radical N | `NG301` | `+0.167` |
| radical O | `OG311` | `-0.379` |
| C2, C6 (quaternary, alpha to N) | `CG301` | `+0.284` |
| C3, C5 (ring CH2) | `CG321` | `-0.237` |
| 4 x gem-dimethyl C | `CG331` | `-0.300` |
| C4 (bears the linker) | `CG314` | `+0.248` |

Bonded terms that must override the release values:

```
BONDS
NG301  OG311   485.17     1.2640   ! Roux
ANGLES
CG301  NG301  CG301    38.60    119.00   ! Roux, near-planar N
CG301  NG301  OG311    45.12    115.52   ! Roux
```

**The residue is +1, and that is not the radical.** The aminoxyl is a neutral
radical: the N-O pair sums to -0.212, a dipole inside a neutral group. The +1
is the protonated secondary amine of the 4-amino-TEMPO linker -- in `RESI LYT`
that nitrogen is `N8`, type `NG3P2`, with **four** bonds (`C23`, `C28`, `H37`,
`H38`), i.e. a dialkylammonium. Radical character and formal charge are
independent properties; TEMPO's own one-electron oxidation product, the
oxoammonium cation, is the charged member of the hydroxylamine / aminoxyl /
oxoammonium triad that this paper's Scheme 1 is built on. The deposited
`TEMPO.pdb` has that nitrogen protonated, so it reproduces the published +1
directly.

`OG311` is borrowed only for its oxygen vdW (`eps -0.1921`, `Rmin/2 1.7650`);
the bonded terms above replace the hydroxyl geometry with nitroxide geometry.
The full set, with provenance and the transfer rule for the charge group, is
in `inputs/nitroxide_overrides.prm`. Load it after
`par_GUI.prm`.

**Citation (confirmed).** D. Sezer, J. H. Freed and B. Roux,
*J. Phys. Chem. B* **2008**, *112*, 5755-5767. Lumata et al. cites this as
ref 71, and its methods read: "employing the CHARMM all-atom force field and
TEMPO parameters from Sezer et al.,⁷¹,⁷² supplemented by the CGenFF force
field.⁷³" Refs 72 and 73 are Vanommeslaeghe/MacKerell, *J. Chem. Inf. Model.*
**2012**, *52*, 3155-3168 and 3144-3154. So the paper already documents the
hybrid Sezer-nitroxide + CGenFF scheme this GUI automates.

This is the general pattern the GUI's CGenFF step is meant to support: where
CGenFF has no chemistry for a group, hand-write the parameters once, keep them
in a documented override file, and let the GUI load them alongside the
machine-generated ones.

### Repairing the input structure

The Open Babel drawing originally supplied as `TEMPO.pdb` had four defects,
all fixed by `inputs/fix_tempo_input.py`, which is idempotent -- it detects and
fixes only what is present, so re-running it is safe. The 2,6-dimethyl original
has been deleted; the script is the record of what was changed.

1. **Two methyls missing.** Each carbon alpha to the nitroxide carried one
   methyl and one hydrogen, making the molecule 2,6-dimethylpiperidin-1-oxyl
   rather than 2,2,6,6-tetramethyl. The gem-dimethyls are what make the radical
   persistent, and they are what make those carbons type as `CG301` so the
   parameters above apply at all. Each alpha C-H was promoted to a methyl
   (C-C 1.530 A, three H at tetrahedral geometry).
2. **Linker amine drawn neutral.** The 4-amino nitrogen had three bonds and one
   hydrogen. The published model has it **protonated** -- a dialkylammonium,
   four bonds, CGenFF type `NG3P2` -- and that is what makes the conjugated
   residue +1. A second hydrogen was added at the fourth tetrahedral vertex and
   the nitrogen marked `N1+`, matching how the 2023 input `tempo_K.pdb` marked
   both of its nitrogens. The structure now reads
   `CC1(C)CC([NH2+]CCCCCC(=O)O)CC(C)(C)N1O`.
3. **Hydroxylamine geometry at N.** The N-O bond was 1.406 A and the nitrogen
   was pyramidal (angle sum 325 deg) -- an sp3 amine, not an aminoxyl. O was
   re-placed in the C-N-C plane at 1.264 A, the override's equilibrium length,
   giving a planar nitrogen (angle sum 360 deg).
4. **Strain from the insertion.** All hydrogens and all four methyl carbons
   were relaxed with MMFF94 while every pre-existing heavy atom (ring, linker,
   acid) was held fixed, so the ring and linker geometry is unchanged.
   Energy 187.3 -> 58.0 kcal/mol; closest non-bonded contact 2.01 A, which is
   the expected 1,3-diaxial methyl H...H of a tetramethylpiperidine.

Result: `C15 H30 N2 O3`, SMILES `CC1(C)CC([NH2+]CCCCCC(=O)O)CC(C)(C)N1O`
(written as the N-OH tautomer only because the sanitiser needs a filled
valence; the deposited file has no hydrogen on that oxygen). That is
4-amino-TEMPO carrying a 6-aminohexanoic acid linker.

The `N1+` formal charge on the nitroxide nitrogen was **deliberately kept**.
The 2023 production input `TEMPO/tempo_K.pdb` carried it too, and that
zwitterion drawing is what the working CGenFF submission used -- a neutral
three-coordinate N with a bare O is a radical, which the CGenFF program will
not accept. The consequence (oxygen typed `OG312`, residue net `-1`) is
corrected by the override file, not by redrawing the input.

**The input now has 49 atoms rather than 43, so CGenFF will return different
atom names. Re-run the CGenFF step and re-take the atom picks above before
using this case study's `*T1`/`*T2`/`*T3` values.**

## Scientific Record

*To be completed by the author:* CGenFF version/date for the `RESI ONE`
parameterisation; the MD protocol used for this case study, if any (engine,
force field files, ensemble, timestep, duration); and whether the solvent-exposure
comparison of Fig. S2 in the publication should be reproduced here and with which
tool.

## Limitations and assumptions

The T158K side chain is modelled, not observed: only `N, CA, C, CB, O` come from
the crystal structure, and `CG, CD, CE, NZ` plus the TEMPO adduct are generated
from internal coordinates. A conformational search was not performed, so the
starting geometry is one plausible rotamer rather than a preferred one.

`RESI ONE` carries CGenFF penalties of 25.0 (parameters) and 34.67 (charges). The
charge penalty is high and concentrated on the nitroxide (`N1` 34.67, `O1`
33.56), which is expected — CGenFF has little coverage of stable radicals. The
nitroxide is the functional centre of the probe, so any quantitative conclusion
about its electronic behaviour should be treated with caution; the parameters are
adequate for packing and solvent exposure, which is what this workflow is for.

The capsid is built as 49 coat-protein copies of one helical segment, not a
complete 300 nm virion, and every surface copy is modified. The publication
reports a distribution of labelling rather than complete occupancy, so the model
represents the fully labelled limit.

## Citation

Lumata, J. L.; Hagge, L. M.; Gaspar, M. A.; Trashi, I.; Ehrman, R. N.; Koirala, S.;
Chiev, A. C.; Wijesundara, Y. H.; **Darwin, C. B.**; Pena, S.; Wen, X.;
Wansapura, J.; Nielsen, S. O.; Kovacs, Z.; Lumata, L. L.; Gassensmith, J. J.
*TEMPO-conjugated tobacco mosaic virus as a magnetic resonance imaging contrast
agent for detection of superoxide production in the inflamed liver.*
J. Mater. Chem. B, 2024, **12**, 3273–3281.
