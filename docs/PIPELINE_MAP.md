# eviVLP chemical-modification pipeline map

File hand-offs for the `chem_GUI` workflow. Working directory for every step is the
folder the GUI runs in (`eviVLP/GUI` locally). Scripts live in `scripts/`, Tcl/psfgen
files in `tcl/`, C++ helpers in `src/` (`IC_naa`, `IC_ntrm_NPC` are compiled from
`IC_table_naa.cpp`, `IC_table_ntrm_NPC.cpp`). Reference PDBs are **not** in git (see end).

Legend: **in →** files read, **→ out** files written. `aa` = `lys` | `cys` | `tyr`.

## The idea in one paragraph

A modified residue is built as one small molecule (model amino acid + drug, with
peptide-bond caps `CX/OX/NX/HX`), parameterized whole by CGenFF, then the caps are cut
off and the protein-side atoms are renamed *by atom name* back to CHARMM36 types.
Parameter lines that become all-CHARMM36 are dropped (redundant), drug-only lines are
kept, and mixed lines are kept with the protein types renamed — those mixed lines are
the cross-force-field parameters. The result is a full `RESI` that *replaces* the
residue in the coat protein, so the model amino acid must be the residue actually
present. For the N-terminus the model is alanine (`BBQ`), which is Qβ's real Ala1.

**Why the protein side is renamed rather than left as CGenFF.** An earlier iteration typed the
entire novel amino acid with CGenFF and handled the cross-force-field terms at the backbone.
That was abandoned for two reasons: (1) each force field's parameters were carefully fitted
as a set, and re-typing protein atoms with CGenFF types would invalidate the CHARMM36
protein parameters (Dr. Nielsen's objection); (2) the N-terminus needed a different scheme
anyway. The current design keeps CHARMM36 types on every protein atom and CGenFF types on
every drug atom, adding only the mixed junction terms. The leftover cross-FF blocks that
used to sit in `par_GUI.prm` came from that earlier iteration; `par_GUI_template.prm` is
the clean base and "Run Renaming" now resets `par_GUI.prm` from it.

## Path A — known residue (Lys / Cys / Tyr)

| # | GUI button / stage | script or program | in → | → out |
|---|---|---|---|---|
| A0 | Load drug | (GUI) copies chosen PDB | drug PDB | `one.pdb`, `one.mol2` (obabel) |
| A1 | CGenFF (web) | user uploads `one.mol2` | `one.mol2` | `one.str` |
| A2 | Orient | `orient_place_drug_aa.py` (filled from `_template`) | `aa.pdb` (capped model residue), `one.pdb` | `one_aligned.pdb` |
| A3 | Patch maker | `patch_maker_aa.py` ← `FileProcessor::ProcessPatchMakerFile` fills `atom_name`/`atom_info` from the 3 bonded + 3 deleted atoms picked in the viewer | `patch_ZAP_aa_template.txt`, `one.str`, `top_all36_cgenff_CBD.rtf` | `PRES ZAP` appended to `top_all36_cgenff_CBD.rtf` |
| A4 | psfgen | `patched_aa.pgn` | `aa.pdb`, `one_aligned.pdb`, rtf | `patched.pdb`, `patched.psf` |
| A5 | element fix | `element_column_fix.py patched.pdb naa.pdb` (resname → NAA) | `patched.pdb` | `naa.pdb` |
| A6 | mol2 | `obabel naa.pdb -O naa.mol2` (no `-h`), then `reconcile_mol2_with_psf.py` (PSF is source of truth for atoms/bonds) | `naa.pdb`, `patched.psf` | `naa.mol2` |
| A7 | CGenFF (web) | user | `naa.mol2` | `naa.str` |
| A8 | Create topology | `noX.py` | `naa.pdb` | `nad.pdb` (cap atoms removed) |
| A9 | | `noX_psf.py` | `naa.pdb`, `naa.str` | `nad.str` |
| A10 | | `rename_atom_types_naa.py aa nad.str` | `nad.str` | `nad.str` (protein types → C36; redundant params dropped; permutations added) |
| A11 | | `fix_atom_naming_pdb_naa.py` | `nad.str`, `nad.pdb` | `nad.pdb` (names from str, resname NAA, resid 1) |
| A12 | | `./IC_naa` | `nad.str`, `nad.pdb`, rtf | `RESI NAA` + IC table appended to `top_all36_cgenff_CBD.rtf` |
| A13 | | `add_naa_params_to_prm.py nad.str <aa> NAA` | `nad.str`, `par_all36m_prot.prm` | novel CGenFF + cross-FF lines into `par_GUI.prm` |
| A14 | Create/view capsid (Qβ) | `FileProcessor::ProcessPgnFile` fills `fix_{a,b,c}_template.pgn` → `fix_{a,b,c}.pgn`; VMD runs them | `qb_cp{1,2,3}_aa.pdb` (coat proteins with target residues pre-renamed `NAA`) | `fix_{a,b,c}.pdb/.psf` |
| A15 | | `Qbeta_capsid_maker.py` → `Qbeta_pgn_writer.py` → `QB_capsid.pgn` | `fix_a.pdb` … | `QB_A_com*.pdb`, `capsid_naa.pdb/.psf` |

NAMD test (`naa_test.conf`) needs `par_all36m_prot.prm`, `par_GUI.prm`, `par_all36_cgenff.prm`.

## Path B — N-terminus (NPC), runs automatically alongside a **Lysine** run

| # | stage | script or program | in → | → out |
|---|---|---|---|---|
| B2 | Orient | `orient_place_drug_ntrm.py` | `BBS_QB.pdb` (Ala1 model, names = `RESI BBQ`), `one.pdb` | `one_aligned.pdb` ⚠ same name as A2's output |
| B3 | Patch maker | `patch_maker_ntrm.py` (filled with the **same** bonded/deleted atoms as the Lys run) | `patch_ZAP_ntrm_template.txt`, `one.str`, rtf | `PRES ZPP` appended to rtf (deletes `1HN2`,`1HN3`; keeps `HN1`; `1CA` = `CG311`) |
| B4 | psfgen | `patched_ntrm.pgn` | `BBS_QB.pdb`, `one_aligned.pdb` | `ntrm.pdb`, `ntrm.psf` |
| B5 | element fix | `element_column_fix_NPC.py ntrm.pdb ntrm_fixed.pdb` (resname → NPC) | `ntrm.pdb` | `ntrm_fixed.pdb` |
| B6 | mol2 | `obabel ntrm_fixed.pdb -O ntrm.mol2` | | `ntrm.mol2` |
| B7 | CGenFF (web) | user | `ntrm.mol2` | `ntrm.str` |
| B8 | Create topology | `noX_ntrm.py` | `ntrm_fixed.pdb` | `ntrm_clean.pdb` (HX removed) |
| B9 | | `noX_psf_ntrm.py` (by **name**; aborts if pdb/str atoms disagree) | `ntrm_fixed.pdb`, `ntrm.str` | `ntrm_clean.str` |
| B10 | | `charge_fix_ntrm.py` (HX charge folded into C) | `ntrm.str`, `ntrm_clean.str` | `ntrm_clean.str` |
| B11 | | `rename_atom_types_naa.py ntrm ntrm_clean.str` (map: N→NH1 CA→CT1 C→C O→O HA→HB1 CB→CT3 HB*→HA3 HN1→H) | `ntrm_clean.str` | `ntrm_clean.str` |
| B12 | | `fix_atom_naming_pdb_ntrm.py` | `ntrm_clean.str`, `ntrm_clean.pdb` | `ntrm_clean.pdb` |
| B13 | | `add_npc_params_to_prm.py` | `ntrm_clean.str` | cross-FF + drug lines into `par_GUI.prm` |
| B14 | | `./IC_ntrm_NPC` | `ntrm_clean.str`, `ntrm_clean.pdb`, rtf | `RESI NPC` + IC table appended to rtf |
| B15 | Capsid | A14 with `{{NTERM}}` → `; first NONE; mutate 1 NPC` (set when `ntrm_clean.str` exists) | `qb_cp*_lys.pdb` | residue 1 built as `NPC` |

Validated end state (May 2026 NAMD run `npc_test`): 36 atoms, net charge 0, single `HN1`
on N, `C` carries the former `HX` charge.

## Junction chemistry (why each file exists)

* `BBS_QB.pdb` / `RESI BBQ`: Ala with protonated amine (`HN1-3`) and an aldehyde cap `HX`
  on `C` so CGenFF sees a closed-shell molecule. `ZPP` removes two amine H's → neutral
  secondary amine N bonded to the drug. `noX` removes `HX`; `charge_fix` keeps the
  residue neutral; `C` then bonds `+N` (Lys2) via `BOND C +N` written by `IC_ntrm_NPC`.
* `CX/OX/NX/HX*` in Path A are the peptide neighbours; same removal logic (`noX*`).
* Patch templates hold 3 bonded slots (`2t*1-3`, drug atoms whose charges are rewritten)
  and 3 deleted slots (`2t*4-6`). **Picking fewer than 3 bonded atoms shifts the deleted
  slots** — a known hazard for new molecules.

## Fixed 2026-10-01/05 (not yet run through VMD/CGenFF)

* `BBS_QB.pdb` rebuilt with `RESI BBQ` names (old Open Babel file kept as `BBS_QB_old.pdb`).
* `patch_ZAP_ntrm_template.txt`: deletes `HN2`/`HN3` (was `HN1`), `1CA` → `CG311`, drug slots templated.
* `noX_psf_ntrm.py`: name-based with consistency check (was index-based; had deleted drug `C1`).
* `fix_atom_naming_pdb_nterm.py` → **`fix_atom_naming_pdb_ntrm.py`** (name the GUI calls), now
  operates on `ntrm_clean.*`.
* `add_npc_params_to_prm.py`: protein-type set matches the `ntrm` map (`CT3`,`HA3`).
* `complete_cross_ff.py` added: copies CGenFF-release values onto renamed protein-side types, closing the junction-parameter gap that was previously patched by hand; GUI now tells the user to enable CGenFF's "include existing parameters" (2026-10-13).
* CGenFF `NBFIX` lines referencing nucleic-acid/carbohydrate/lipid types stripped by `filter_cgenff_nbfix.py` (NAMD cannot skip them the way CHARMM's BOMLEV can) (2026-10-13).
* **Root cause of the NAMD failures: a CGenFF v5.0 `par_all36_cgenff.prm` installed 2026-05-15 that NAMD 3.0b5 cannot read, beside a v4.6 topology. Restored v4.6.**
* `par_GUI.prm`: vestigial CMAP section (renamed CGenFF-typed copy of the CHARMM36 backbone correction) deleted - six extra maps on top of the standard six made NAMD fail; orphaned-CMAP strip bug repaired (a bad strip had eaten the `CMAP` header) and a terminating `END` added; `validate_prm.py` now gates it (2026-10-13).
* `par_GUI.prm` given a terminating `END` (NAMD died with ABNORMAL EOF whenever it was the last parameter file); parameter scripts now stop inserting at NONBONDED/END (2026-10-13).
* IC generators now walk bonded paths only (1207 -> 266 entries for RESI NAA); junction-angle distortion in the orient scripts diagnosed and logged (2026-10-13).
* NAMD log redirect fixed (`> log 2>&1`; the log had been passed as a second config-file argument, so NAMD wrote nothing); failures now print the log tail (2026-10-13).
* Capsid stage split into Build Monomer / Minimize Monomer / Assemble Capsid, since the NAA residue is constructed during the monomer build and must exist before it can be minimised (2026-10-13).
* Monomer minimisation stage added (`write_fixed_atoms.py`, `minimize_monomer_template.conf`, `coor_to_pdb_template.tcl`, "Minimize Monomer" button) (2026-10-08).
* A13 wired: `add_naa_params_to_prm.py` now called by the GUI, reading protein types from the rename map and checking CHARMM36 before dropping a line; output zip follows the VLP selection (2026-10-08).
* Charge overrides added: `apply_charge_overrides.py` plus an "Apply charge overrides" button (optional, beside the NAA CGenFF import). CGenFF types a persistent radical's oxygen as an alkoxide (`OG312`), which also gives the residue a formal charge it should not have, so types and charges are replaced from a data file after "Run Renaming" and before psfgen. The script refuses when the residual per balance atom exceeds 0.25 e, since that means either the override under-specifies the charge group or the molecule's protonation state does not match the target - neither is fixable by redistributing charge. Nitroxide citation confirmed: Sezer, Freed & Roux, J. Phys. Chem. B 2008, 112, 5755-5767 (ref 71 of Lumata et al. 2024) (2026-10-13).
* Base GUI no longer ships a case-specific construct: lysine on TMV now resolves to the native `tmv_K53_K68.pdb`; `tmv_K158.pdb` (T158K mutant) moved to the case study inputs. `.gitignore` excepted `case_studies/**/inputs/*.pdb` but case study 1 used `input/`, so none of its PDBs had ever been committed - folder renamed to `inputs/` (2026-10-13).
* TEMPO input structure repaired: the drawing was 2,6-dimethyl, not 2,2,6,6-tetramethyl, and had hydroxylamine geometry at N; `fix_tempo_input.py` added. Aminoxyl radical parameters recovered from the 2023 runs and documented in `nitroxide_overrides.prm` - CGenFF has no radical type in v4.6 or v5.0 (2026-10-13).
* Case study 1 (kTMV/T158K) documented; `tmv_K158.pdb` added (2026-10-08).
* TMV templated for lysine and a Ser N-terminus model (`BBS_TMV.pdb`/`BBT`, `ntrm_ser`) added; `ntrm_ref.pdb` indirection (2026-10-08).
* IC generators de-duplicated; `validate_rtf_ic.py` added; amide-H IC name corrected (2026-10-08).
* CGenFF import verifies the .str against its mol2 (`check_str_matches_mol2.py`); residue name set in the mol2; Patch Maker / Create NAA re-runnable via `strip_rtf_residues.py`; VMD success judged by output file (2026-10-08).
* `rename_atom_types_naa.py`: amide H type `H` (was invalid `HN`); acylated NZ -> `NH1`; all-protein parameter lines only dropped when `par_all36m_prot.prm` really defines them (2026-10-08).
* `par_GUI_template.prm` created; `par_GUI.prm` reset from it at "Run Renaming" (2026-10-08).
* Orient step: auto-run from Patch Maker when `one_aligned.pdb` is stale; orient scripts exit 1 on error; N-terminus orient script now gets its placeholders filled (2026-10-08).
* CGenFF import accepts the result zip; mol2 download via save dialog; startup python3/obabel check (2026-10-08).
* `fix_{a,b,c}_template.pgn` created; `FileProcessor::ProcessPgnFile` fills all three and the
  `{{NTERM}}` slot; `chem_GUI.cpp` now checks its return value. **GUI must be recompiled.**

## The degenerate-IC bug (fixed 2026-10-08)

`writeInternalCoordinateTable` in `IC_table_naa.cpp` / `IC_table_ntrm_NPC.cpp` looped over
four independent indices with no requirement that they be distinct, so ~53% of the IC
entries it wrote had a repeated atom (`IC C1 O1 C1 H16`, `IC N H N CA`, angle ~0).
psfgen's `guesscoord` uses the FIRST IC it can satisfy for an atom; a degenerate entry
needs fewer distinct known atoms than the correct one, so it won and placed the atom
from a zero-degree frame. Symptom: a drug hydrogen (e.g. the ring N-H) lands on top of
the backbone amide H, 13 A from the nitrogen it is actually bonded to, so in VMD it
looks as though a backbone hydrogen has swapped its bond to a ring nitrogen. The PSF
topology was always correct; only the coordinates were wrong.

Fixes: a distinctness guard in both generators (1207 -> 564 IC lines for RESI NAA, 0
degenerate); the hard-coded trailer IC now names the amide hydrogen `H` for NAA and
`HN1` for NPC (it said `HN`, a CHARMM22/27 name absent from both residues); and
`validate_rtf_ic.py` runs after each IC step, failing the GUI step on a degenerate
entry, an undefined atom name, or an atom with no usable IC.

## VLP coverage and the two N-terminus models

Both VLPs now use the same template mechanism:

* **Qbeta** - `fix_{a,b,c}_template.pgn`, slots `{{AMINO_ACID}}` and `{{NTERM}}`.
  Coat proteins `qb_cp{1,2,3}_{lys,cys}.pdb`. Residue 1 is **Ala**.
* **TMV** - `TMV_naa_template.pgn`, slots `{{COAT_PDB}}` and `{{NTERM}}`, filled by
  `FileProcessor::ProcessPgnFile` -> `TMV_naa.pgn`. Coat proteins are resolved by
  `GetTmvCoatProtein()`: tyr -> `tmv_Y139.pdb` (Tyr139, nTMV diazonium + CuAAC),
  lys -> `tmv_K53_K68.pdb` (the **native** Lys53/Lys68 doxorubicin sites of
  Finbloom et al. 2016; both residues are marked `NAA`, so a run conjugates at
  both). No cysteine template exists; the GUI says so and stops.

  The native sites keep their crystallographic `CG/CD/CE/NZ`, so `guesscoord`
  builds only the drug, and residue 158 stays `THR` so `CTER` is applied to an
  ordinary residue. The psfgen templates are residue-count agnostic - they load
  the whole coat protein as segment `A` and build every `NAA` from the same
  topology - so one site and two behave identically.

  **Case study 1 (kTMV) is deliberately not the default.** `tmv_K158.pdb` now
  lives in `case_studies/case_study1_TMV_TEMPO/inputs/`: it is the T158K
  single-point mutant of Lumata et al. JMCB 2024, so wild-type TMV has no lysine
  there, and it keeps only `N, CA, C, CB, O` from Thr158 - the whole lysine side
  chain plus the adduct come from `guesscoord`, making the conformation a model
  rather than an observation. **Residue 158 is also the C-terminus**, so psfgen
  applies its default `CTER` patch to the `NAA` residue itself: `O` deleted,
  `OT1`/`OT2` added, `C` retyped `CC`. `RESI NAA` supplies the `CA`, `C`, `O` that
  `CTER` needs, and the dangling `BOND C +N` from `IC_naa` is ignored. That is the
  first terminal NAA the toolkit has built and it has not been through psfgen yet.
  To reproduce that arm, put the file beside the GUI and return its name from
  `GetTmvCoatProtein()`.
  `TMV_pgn_writer.py` adds the same `; first NONE; mutate 1 NPC` per segment when
  `ntrm_clean.str` is present. Residue 1 is **Ser**.

**Why the N-terminus is modified at all.** Amine-directed chemistry (NHS esters,
activated acids) cannot distinguish the lysine side-chain epsilon-amine from a free
N-terminal alpha-amine, so a lysine job necessarily modifies both - which is why the
NPC path runs alongside every Lysine run. This assumes the alpha-amine is FREE.
Recombinant coat proteins expressed in *E. coli* are not N-terminally acetylated
(*E. coli* lacks the eukaryotic NatA/NatB/NatC acetyltransferases), so for the
E. coli-expressed VLPs this toolkit targets, residue 1 presents a free alpha-amine.
A capsid from a native plant or eukaryotic source may carry an acetylated N-terminus
and be unreactive; the N-terminus should then not be modified at all. The GUI states
this assumption in the status log when it stages the model residue.

The N-terminus model residue must be the amino acid actually at position 1, because
the NPC path replaces that residue whole (see the cross-FF paragraph above):

| VLP | residue 1 | model PDB | RESI | rename map |
|---|---|---|---|---|
| Qbeta | Ala | `BBS_QB.pdb` | `BBQ` | `ntrm` |
| TMV | Ser | `BBS_TMV.pdb` | `BBT` | `ntrm_ser` |

`BBS_TMV.pdb` is derived from `BBS_QB.pdb`'s validated backbone with `HB3` replaced by
`OG` (1.42 A) plus `HG1` (0.96 A, 108 deg); types follow `RESI SER` in
`top_all36_prot.rtf`. All scripts read a single staged filename `ntrm_ref.pdb`;
`MyFrame::StageNtrmReference()` copies the right model there from the VLP selection,
and `NtrmVariant()` supplies `ntrm` / `ntrm_ser` to `rename_atom_types_naa.py` and
`add_npc_params_to_prm.py` (whose protein-type set is now read from the rename map
rather than hard-coded, so the two cannot drift).

**One-time setup per model residue:** `RESI BBT` must exist in
`top_all36_cgenff_CBD_template.rtf`. Send `BBS_TMV.mol2` to CGenFF once, then
`python3 install_nterm_reference.py BBT.str` inserts the topology block (psfgen needs
topology only; the conjugate's parameters still come from the NPC stream file).

Adding a third VLP whose residue 1 is something else means one more model PDB, one more
map entry, and one more CGenFF run - which is why the PRES first-patch design stays on
the list as the general answer.

## What `par_GUI.prm` must contain

A NAMD run loads three parameter files: `par_all36m_prot.prm` (CHARMM36 protein),
`par_all36_cgenff.prm` (released CGenFF), and `par_GUI.prm`. The first two are fixed, so
`par_GUI.prm` carries exactly the remainder, and both paths write into that one file:

* **novel CGenFF terms** for the drug. A CGenFF stream file's parameter section only ever
  holds terms absent from the released set, so every drug-only line qualifies.
* **cross-force-field terms** at the junction — lines mixing CHARMM36 protein types with
  CGenFF drug types, which neither force field defines.
* **CHARMM36 gaps** — an all-protein line is dropped as redundant only when
  `par_all36m_prot.prm` actually defines it, checked directly rather than assumed.

Protein-side types come from the rename map in `rename_atom_types_naa.py`, so the
parameter scripts and the renaming step cannot drift apart. Both scripts replace their own
previous block, so the steps are re-runnable. The output zip ships `par_GUI.prm` and
`top_all36_cgenff_CBD.rtf` together with the solvated capsid.

## Minimisation: two stages, different scopes

The capsid stage is split into three buttons, because the coat-protein templates only
LABEL the target residues `NAA`/`NPC` - the real residue, with its adduct, is constructed
by psfgen from `RESI NAA`/`RESI NPC` plus `guesscoord`. That construction has to finish
before there is anything to minimise, and the minimisation has to finish before the
symmetry matrices replicate the subunit:

| Button | Does | Writes |
|---|---|---|
| 1. Build Modified Monomer(s) | `ProcessPgnFile` then psfgen per monomer | `fix_{a,b,c}.pdb/.psf` (Qbeta) or `TMV_naa.pdb/.psf` |
| 2. Minimize Monomer (protein fixed) | stage 1 below | `<stem>.pdb` updated, `<stem>_premin.pdb` kept |
| 3. Assemble And View Capsid | replicate, assemble, solvate, view | `capsid_naa.*` / `TMV_naa_capsid.*` |

Button 3 builds the monomers only if they are **absent**; an existing monomer may have
been minimised, so it is never rebuilt over. It reports how many of the subunits carry a
`_premin.pdb`, i.e. have been minimised. Conversely button 1 deletes any stale
`_premin.pdb`/`_min.coor`, because a fresh build invalidates an earlier minimisation.

**Stage 1 - monomer, protein fixed, BEFORE replication** ("Minimize Monomer"):
`write_fixed_atoms.py <stem>.pdb <stem>_fix.pdb --only-residues NAA NPC` →
`minimize_monomer_template.conf` → NAMD (`minimize 1000`, non-periodic, no solvent) →
`coor_to_pdb_template.tcl` writes the result back over `<stem>.pdb`
(original kept as `<stem>_premin.pdb`). Stems: `fix_a`/`fix_b`/`fix_c` for Qbeta,
`TMV_naa` for TMV.

Which atoms move is decided by psfgen's own bookkeeping: `writepdb` sets occupancy
0.00 for coordinates it GUESSED and 1.00 for coordinates from the input structure. So
the free set is exactly what the pipeline invented - the adduct, side-chain atoms a
mutation added, and placed hydrogens - and every crystallographic heavy atom is held.
The script refuses to run on a file whose occupancies are all 1.00 rather than guess.

This matters because the capsid is assembled by applying crystallographic symmetry
matrices to the monomer, and those matrices are valid only for the conformation they
were derived from. If crystallographic atoms relax here, every subunit interface in the
assembled capsid is placed wrong - and an isolated subunit drifts AWAY from its
in-capsid conformation, which is held by inter-subunit contacts that do not exist for a
lone monomer. It matters most for kTMV/T158K, where only `N, CA, C, CB, O` are
crystallographic and the whole side chain plus adduct come from the IC table.

**Stage 2 - assembled capsid, AFTER replication** ("Minimize Capsid", `minimize.conf`):
the solvated, ionised capsid, everything free. Only this stage can reveal clashes
BETWEEN neighbouring subunits' adducts, which is a real possibility because the adduct
sits on the solvent-exposed surface. Running it on raw `guesscoord` geometry means
hundreds of bad contacts at once, which minimisers handle badly and which can distort
protein that was fine - hence stage 1 first.

`minimize.conf` is Qbeta-specific (`QB_naa_capsid_wb`, 320 A cell, `cellOrigin 74 0 0`)
and the GUI now warns before running it on a TMV selection. A TMV run needs the
structure/coordinates names and the cell set for the rod.

Either stage doubles as a parameter check: a missing parameter aborts NAMD during setup,
before step one, so stage 1 answers "is `par_GUI.prm` complete?" in seconds without
building a capsid.

## IC frames are built from bonded paths (fixed 2026-10-13)

`writeInternalCoordinateTable` accepts any quadruple A-B-C-D whose consecutive distances
are under 2.1 A. That threshold admits pairs that are close but **not bonded**: two
hydrogens on the same carbon sit 1.79 A apart, and a geminal H...C is about 2.08 A. The
result is IC entries built on degenerate frames, e.g.

    IC CA   CB   HB1  HB2   1.5573 108.93  -117.845 36.0385  1.79367
    IC CG   CB   HB2  HB1   1.5369 106.643 -116.835 36.1961  1.79367

a 36-degree angle with an H-H "bond length". psfgen then reports
`failed to guess coordinate due to bad angle CG CB HB2` and
`poorly guessed coordinate for atom HB2`. The distinctness guard added earlier removed
entries with a *repeated* atom; these have four distinct atoms in a bad geometric
arrangement, so they survive it.

**Fixed:** `buildAdjacency()` parses the RESI block's `BOND`/`DOUBLE`/`TRIPLE` lines into
an adjacency list, and the generator now walks genuine bonded paths A-B-C-D instead of
every quadruple. Frames with a near-collinear A-B-C or B-C-D (under 5 or over 175 degrees)
are also skipped. Note `readData` keeps the raw 4-character PDB name field, so
`atoms[].type` carries padding and name lookups go through `trimName()`; the IC output
keeps the padded form so the table stays column-aligned.

Measured on `RESI NAA` (61 atoms, 60 bonds): **1207 -> 564 -> 266 IC entries**
(original, after the distinctness guard, after the bonded-path requirement), zero
degenerate frames, zero H-H frames, and the validator reports every atom placeable. `HB2`
is now placed from `IC N CA CB HB2 ... 109.2 ... 108.8` instead of a 36-degree H-H frame.

Generation cost drops from O(N^4) to O(N x degree^3), which matters for a larger adduct:
61 atoms was 13.8M quadruple evaluations, now a few hundred paths.

## The amide junction geometry is distorted (open)

`orient_place_drug_<aa>.py` aligns the drug by making the `*T1 -> *T3` vector parallel to
the residue's `CD -> NZ` vector and placing `*T1` at the target distance along that same
line. That fixes the bond length correctly but leaves `CD`, `NZ`, `*T1` and `*T3`
**collinear**, which is wrong for an amide: both the nitrogen and the carbonyl carbon
should be roughly trigonal.

Measured on the last lysine run's `nad.pdb`, bond lengths are right and angles are not:

| | measured | ideal |
|---|---|---|
| `NZ-C12` | 1.350 A | 1.35 |
| `C12=O6` | 1.218 A | 1.23 |
| `O6-C12-NZ` | **51.5 deg** | ~123 |
| `CE-NZ-C12` | **144.5 deg** | ~120 |
| `NZ-C12-C11` | **144.8 deg** | ~116 |

The three angles around `C12` sum to 322 degrees rather than ~360, so the carbonyl carbon
is pyramidalised with `O6` folded back toward the nitrogen. This is the same underlying
cause as the earlier "the oxygen is closer than the carbon" problem: `reconcile_mol2_with_psf.py`
fixed the *bonds*, but the *geometry* is still wrong.

The monomer minimisation stage is the current mitigation and is well motivated by this -
the bonds and types are correct, so a minimiser can recover the angles. A proper fix is a
real geometric construction in the orient script: place `*T1` along an sp3 direction from
`NZ` and then rotate the adduct so `*T1`'s substituents splay at ~120 degrees about the
`NZ-*T1` bond, rather than aligning `*T1 -> *T3` to the `CD -> NZ` axis.

## Normal psfgen messages (not errors)

`skipping conformation C-N-CA-C at beginning of segment`, `skipping bond C-N at end of
segment`, `skipping improper C-CA-N-O at end of segment`, `skipping cross-term ... at end
of segment` are all expected: they are the IC/BOND/IMPR entries that reference the
previous (`-C`) or next (`+N`) residue, which do not exist at a chain terminus. psfgen
correctly drops them.

`failed to set coordinate for atom O TYR:132` means that backbone oxygen is **absent from
the input coat-protein PDB**, so psfgen guesses it. Worth knowing because that atom is
modelled rather than observed; it is not caused by this pipeline.

## Enable CGenFF's "include parameters already in CGenFF" (2026-10-13)

The CGenFF web service offers a checkbox for whether the returned stream file should
include parameters it **already had**, as well as the new analogy-derived ones. It is OFF
by default, and the files in the GUI folder were generated that way: `naa.str` and
`ntrm.str` each carry just **2** parameter lines, both with nonzero penalties.

With it off, the stream file is not self-contained. Every term the molecule needs beyond
those two must come from the locally installed `par_all36_cgenff.prm`, so the run depends
on the local release matching the service's library. When it does not, terms have no
source anywhere - which is the failure that led to a CGenFF v5.0 file being installed
beside a v4.6 topology in May, and everything that followed from that.

**Turn it on.** Then the stream file carries every parameter the molecule needs, with the
values the service actually assigned, and `complete_cross_ff.py` can take donors from it
rather than from the local release. The script prefers the stream file for exactly that
reason and reports donor provenance; it also prints a NOTE when a stream file looks like
it was generated with the option off. The costs are benign: more `DUPLICATE DIHEDRAL`
warnings from NAMD and a larger `par_GUI.prm`, and since `par_GUI.prm` is read last its
values win.

## Renaming orphans terms CGenFF took from its own library (fixed 2026-10-13)

`UNABLE TO FIND ANGLE PARAMETERS FOR CG321 CG2O1 NH1`. CGenFF's release defines
`CG321 CG2O1 NG2S1` (80.00, 116.50); the rename maps the backbone N `NG2S1 -> NH1`, so
the structure needs `CG321 CG2O1 NH1`, which exists nowhere. A stream file only contains
terms CGenFF could **not** find in its own library, so a term it *did* find never appears
there and nothing re-emits it under the renamed type. The `**cgenff_prot_backbone` angles
already in `par_GUI.prm` were this class patched by hand - `CG311 CG2O1 NH1` was present,
the `CG321` variant was not.

`complete_cross_ff.py` closes it generally: enumerate every bonded term the residue
contains (bonds, angles, dihedrals, impropers from the BOND/IMPR lists), check whether the
**renamed** type combination is available in the loaded set, and where it is not, copy the
numbers from the equivalent **original-type** term. It only copies values a force field
already defines and never invents one. Measured: NPC 1 angle + 5 dihedrals + 1 improper;
NAA 15 angles + 6 dihedrals + 1 improper - so this was a systematic gap, not a one-off.

Two traps found while building it: generated lines need trailing newlines or several
dihedrals end up on one line; and the availability check must exclude the block it is
about to replace, or a re-run fills only the newly-missing terms and then wipes the rest.
Both are guarded, and two consecutive passes now produce an identical file.

## CGenFF NBFIX references types from force fields we do not load (fixed 2026-10-13)

After the version fix below, NAMD got through reading and then died with

    FATAL ERROR: DIDN'T FIND vdW PARAMETER FOR ATOM TYPE OC2D1

`OC2D1` appears in **no** atom of the structure and in **no** topology file - it is
referenced only by CGenFF's own `NBFIX` section, which carries pair corrections against
the nucleic-acid, carbohydrate and lipid force fields. 17 such types are unsatisfied when
only protein + CGenFF + par_GUI are loaded: `NN1 NN2 NN2G NN2U NN3 NN3A NN3G NN4 ON1
ON1C` (nucleic acid), `OC2D1-OC2D4` (carbohydrate), `OBL OCL NH3L` (lipid). The CGenFF
header documents the condition and tells CHARMM users to set `BOMLEV (-2)`; NAMD has no
equivalent and treats it as fatal.

A VLP capsid in water with ions contains no nucleic acid, carbohydrate or lipid, so those
corrections can never apply. `filter_cgenff_nbfix.py` writes a copy of the released file
with exactly those `NBFIX` lines removed (33 of them) and a banner recording the command
that produced it:

    python3 filter_cgenff_nbfix.py par_all36_cgenff_release_v4.6.prm \
            par_all36_cgenff.prm par_all36m_prot.prm par_GUI.prm

The released file is kept untouched as `par_all36_cgenff_release_v4.6.prm`, so the filtered
one is reproducible and provenance is preserved. This is the one place the vendor file is
derived rather than used as-is, and it is a deletion of inapplicable terms - never an
addition. Novel and cross-force-field parameters still belong in `par_GUI.prm`.

**A CGenFF stream file never contains a NONBONDED section** (verified: zero in `one.str`,
`nad.str`, `ntrm_clean.str`). It supplies bonded terms only. So every atom type the CGenFF
service assigns must already have a vdW entry in the installed parameter file - which is
the hard reason the topology and parameter versions have to match, and why a new type
cannot be rescued by the stream file alone.

## The ABNORMAL EOF was par_all36_cgenff.prm, not par_GUI.prm (fixed 2026-10-13)

Isolating the parameter files one at a time settled it:

| files loaded | result |
|---|---|
| prot + cgenff | `ABNORMAL EOF FOUND - buffer=*END*` |
| prot only | reads fine, then the expected `DIDN'T FIND vdW PARAMETER FOR ATOM TYPE CG331` |
| **cgenff only** | `ABNORMAL EOF` |

`par_all36_cgenff.prm` in the GUI folder was a **CGenFF v5.0** file (3.6 MB, 29041
dihedrals) that NAMD 3.0b5 cannot read, and it was also a version mismatch: the topology
`top_all36_cgenff_CBD_template.rtf` is **v4.6**, and the stream files the CGenFF web
service returns say *program version 4.0*. File timestamps date the regression exactly:

| | |
|---|---|
| `npc_test.log` - last successful run | 2026-05-14 14:33 |
| `par_all36_cgenff.prm` v5.0 installed | 2026-05-15 16:08 |
| `naa_test.log` - first failure | 2026-05-15 16:08 |

Every run since has failed the same way; the error only looked new because nothing had
reached NAMD in months. Restored **v4.6** from
`junk_maybe/silcsbio.2024.1/data/cgenff/par_all36_cgenff.prm` (the copy that also carries
the SOD/CLA ion parameters); the v5.0 file is kept as `par_all36_cgenff_v5.0.prm`. All 17
atom types used by `nad.str` and `ntrm_clean.str` have vdW parameters in 4.6 - no gaps.

**Keep the released CGenFF file unmodified.** Novel and cross-force-field parameters
belong in `par_GUI.prm`, read **last** so they take precedence; a later file overrides an
earlier one for the same term. The `DUPLICATE DIHEDRAL ENTRY` / `INCREASING MULTIPLICITY`
warnings are expected in that arrangement and are how our values win. Editing the vendor
file instead would lose provenance and make a run impossible to reproduce.

## par_GUI.prm: two further defects found along the way (both fixed 2026-10-13)

### 0. The CMAP section should not be there at all (the actual cause)

`par_GUI.prm` carried six **CMAP** maps that are a renamed copy of the CHARMM36 backbone
correction, with CGenFF types at the termini:

    C  NH1  CT1  C  NH1  CT1  C  NG2S1   24
    CG2O1  NH1  CT1  C  NH1  CT1  C  NH1   24

They are a leftover from the earlier all-CGenFF novel-amino-acid design. Under the
current design the NAA/NPC backbone keeps **standard** CHARMM36 types - `nad.str` has
`N NH1`, `CA CT1`, `C C`, `O O`, and `CG2O1`/`NG2S1` never appear as backbone types - so
the CMAP already in `par_all36m_prot.prm` applies and these six are redundant.

Declaring six extra maps on top of the standard six breaks NAMD: it then dies with
`ABNORMAL EOF FOUND - buffer=*END*` at the end of whichever parameter file is read
**after** `par_GUI.prm`, or at the end of `par_GUI.prm` itself when that is last. That is
why reordering moved the error rather than removing it.

This also explains why runs worked in May: the `CMAP` header was missing then (see
defect 1 below), so NAMD never parsed the block - it absorbed those ~1030 lines as
garbage improper entries instead. Restoring the header, intended as a repair, is what
made NAMD actually read the maps and fail. The section is now **deleted** from
`par_GUI_template.prm`, with a comment in its place recording why.

Parameter order: `par_GUI.prm` is read **last**, because a later file overrides an earlier
one for the same term and our novel/cross-force-field values must win. The log showed
cgenff replacing them when read afterwards (`IDENTICAL PERIODICITY! REPLACING OLD VALUES`).

### 1. An orphaned CMAP block (a regression, now moot)

`par_GUI.prm` carries a **CMAP** cross-term block - 8 atom types followed by a 24x24 grid
of bare numbers. The first attempt at building `par_GUI_template.prm` stripped each
`! === ... cross-force-field ... ===` block by consuming lines until the next **blank**
line, and the `IMPROPERS` block ran straight into the `CMAP` header with no blank line
between. So the header was swallowed and ~1030 lines of grid data were left orphaned
inside `IMPROPERS`. NAMD reads that as malformed improper entries and runs to EOF:

    FATAL ERROR: ABNORMAL EOF FOUND - buffer=*END*

The stripper now consumes only parameter lines and **stops at** a blank line, a section
keyword, or another block header without consuming it. The CMAP section it had damaged is
in any case deleted now (above), so the structure is
`BONDS, ANGLES, DIHEDRALS, IMPROPERS, NONBONDED, END`.

### 2. No terminating END (pre-existing)

Separately, `par_GUI.prm` has never had a terminating `END`, unlike `par_all36m_prot.prm`,
`par_all36_cgenff.prm` and `par_water_ions.inp`. NAMD tolerates that only while another
parameter file is read afterwards; read **last**, it hits EOF still expecting data and
dies with the same message. The two May test logs differ only in that order:

| log | parameter order | outcome |
|---|---|---|
| `naa_test.log` | cgenff, prot, **par_GUI last** | `ABNORMAL EOF` |
| `npc_test.log` | prot, **par_GUI middle**, cgenff | ran: 500 min + 1000 steps |

Note what this means for `npc_test`: it survived because the next file ended the parse,
but the orphaned CMAP lines were still being absorbed as garbage improper parameters. Any
result from a run that loaded this file is suspect.

### Guard

`validate_prm.py` checks section structure, per-section field counts (CMAP grid data
recognised, so it is not mistaken for malformed impropers) and the terminating `END`. The
GUI runs it after inserting parameters and refuses to continue if it fails, so a bad file
can no longer reach NAMD.

### Related latent bug

`par_GUI.prm` had **no terminating `END`**, unlike `par_all36m_prot.prm`,
`par_all36_cgenff.prm` and `par_water_ions.inp`, which all have one. NAMD tolerates that
only while another parameter file is read afterwards; if `par_GUI.prm` is the LAST
`parameters` line, NAMD reaches EOF still expecting data and dies with

    FATAL ERROR: ABNORMAL EOF FOUND - buffer=*END*

The two May test logs in the GUI folder are the proof, and they disagree only in order:

| log | parameter order | outcome |
|---|---|---|
| `naa_test.log` | cgenff, prot, **par_GUI last** | `ABNORMAL EOF` |
| `npc_test.log` | prot, **par_GUI middle**, cgenff | ran: 500 min + 1000 steps |

So this was latent from the start and surfaced as soon as a config listed `par_GUI.prm`
last - which `minimize_monomer_template.conf` did. `END` is now appended to
`par_GUI_template.prm` (the source of truth, copied over `par_GUI.prm` at *Run Renaming*),
so parameter order no longer matters.

Both parameter scripts computed a section's end
using only `CMAP` as a stop, and `par_GUI.prm` has no `CMAP`. An `IMPROPERS` block would
therefore have been appended at end-of-file - after `NONBONDED`, and now after `END`. Both
scripts now stop at the first of `CMAP`, `NONBONDED`, `NBFIX`, `HBOND`, `END`. Only
`DIHEDRALS` blocks have appeared so far, which is why this never bit.

## Still open

* Re-run CGenFF for the current drug: `ntrm.str` in `GUI/` is from the May molecule, `one.pdb` from June.
* `Qbeta_pgn_writer.py` writes `first NONE` for every segment; check whether residue-1 NTER
  hydrogens from `fix_a.pdb` survive the capsid build.
* Generalising the N-terminus beyond Ala1 would mean a `PRES`-based first-patch design
  (the March `patch_NPC_ntrm.txt` idea) with CB-type variants — future work.

## Files the GUI needs that git does not track

`BBS_QB.pdb`, `lys.pdb`, `cys.pdb`, `tyr.pdb`, `qb_cp{1,2,3}_{lys,cys,tyr}.pdb`,
`patch_ZAP_{lys,cys,tyr,ntrm}_template.txt`, `reconcile_mol2_with_psf.py`,
`top_all36_prot.rtf`, `par_all36m_prot.prm`, `par_all36_cgenff.prm`.
`.gitignore` already allows PDBs under `chemical_modification/data/`.
