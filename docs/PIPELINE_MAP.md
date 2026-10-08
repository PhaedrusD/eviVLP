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
| A13 | *(not wired)* | `add_naa_params_to_prm.py` | `nad.str` | cross-FF lines into `par_GUI.prm` — **currently a manual step** |
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
* TMV templated for lysine (K53/K68) and a Ser N-terminus model (`BBS_TMV.pdb`/`BBT`, `ntrm_ser`) added; `ntrm_ref.pdb` indirection (2026-10-08).
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
  `GetTmvCoatProtein()`: tyr -> `tmv_Y139.pdb` (Tyr139), lys -> `tmv_K53_K68.pdb`
  (Lys53 + Lys68). No cysteine template exists; the GUI says so and stops.
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

## Still open

* Re-run CGenFF for the current drug: `ntrm.str` in `GUI/` is from the May molecule, `one.pdb` from June.
* Path A cross-FF parameters are not inserted into `par_GUI.prm` by the GUI (A13).
* `Qbeta_pgn_writer.py` writes `first NONE` for every segment; check whether residue-1 NTER
  hydrogens from `fix_a.pdb` survive the capsid build.
* Generalising the N-terminus beyond Ala1 would mean a `PRES`-based first-patch design
  (the March `patch_NPC_ntrm.txt` idea) with CB-type variants — future work.

## Files the GUI needs that git does not track

`BBS_QB.pdb`, `lys.pdb`, `cys.pdb`, `tyr.pdb`, `qb_cp{1,2,3}_{lys,cys,tyr}.pdb`,
`patch_ZAP_{lys,cys,tyr,ntrm}_template.txt`, `reconcile_mol2_with_psf.py`,
`top_all36_prot.rtf`, `par_all36m_prot.prm`, `par_all36_cgenff.prm`.
`.gitignore` already allows PDBs under `chemical_modification/data/`.
