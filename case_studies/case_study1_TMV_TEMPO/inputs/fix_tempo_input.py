"""Repair the eviVLP TEMPO input structure.

Idempotent: each defect is detected and fixed only if present, so re-running on
an already-corrected file is safe and changes nothing but the MMFF relaxation.

Four defects in the Open Babel drawing, in order of consequence:
  1. Both carbons alpha to the nitroxide carry one methyl and one H, so the
     molecule is 2,6-dimethylpiperidin-1-oxyl, not 2,2,6,6-tetramethyl (TEMPO).
     The gem-dimethyls are what make the radical persistent, and they are what
     make the alpha carbons type as CG301 so the hand-curated nitroxide
     parameters apply. Each alpha C-H is promoted to a methyl.
  2. The 4-amino linker nitrogen is drawn NEUTRAL (three bonds, one H). In the
     published model it is PROTONATED - a dialkylammonium, four bonds, type
     NG3P2 - which is what makes the conjugated residue +1. A second H is added
     and the nitrogen marked N1+ so CGenFF types it as an ammonium.
  3. The nitroxide N is pyramidal (angle sum 325 deg) and the N-O bond is
     1.406 A, i.e. a hydroxylamine. An aminoxyl N is near-planar with
     N-O ~1.26 A. O is re-placed in the C-N-C plane at the override's
     equilibrium length.
  4. Geometry around the new atoms is strained after insertion, so the
     hydrogens and all four methyl carbons are relaxed with MMFF94 while every
     pre-existing heavy atom (ring, linker, acid) is held fixed.

WHERE THE +1 COMES FROM, and where it does not: the aminoxyl is a NEUTRAL
radical - its N-O pair sums to -0.2120, a dipole inside a neutral group. The
residue's +1 is defect 2 above, the protonated linker amine. Radical character
and formal charge are independent properties; TEMPO's own one-electron
oxidation product, the oxoammonium cation, is the charged member of the
hydroxylamine / aminoxyl / oxoammonium triad. Charges follow Sezer, Freed &
Roux, J. Phys. Chem. B 2008, 112, 5755-5767, as used by Lumata et al.,
J. Mater. Chem. B 2024, 12, 3273-3281 (ref 71).

The N1+ formal charge on the NITROXIDE nitrogen is deliberately KEPT: the 2023
production input (TEMPO/tempo_K.pdb) carried it on both nitrogens, and that
zwitterion drawing is what the working CGenFF submission used. The resulting
OG312 typing and net -1 are corrected afterwards by nitroxide_overrides.prm
(alongside this script), not here.

usage: python3 fix_tempo_input.py in.pdb out.pdb
"""
import sys, math, collections
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

CC, CH, NO, NH = 1.530, 1.090, 1.2640, 1.0100
COS_TET, SIN_TET = -1.0/3.0, math.sqrt(8.0)/3.0
sub=lambda a,b:[x-y for x,y in zip(a,b)]; add=lambda a,b:[x+y for x,y in zip(a,b)]
mul=lambda a,s:[x*s for x in a]; dot=lambda a,b:sum(x*y for x,y in zip(a,b))
nrm=lambda a: math.sqrt(dot(a,a)); unit=lambda a: mul(a,1.0/nrm(a))
def cross(a,b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

src,dst = sys.argv[1], sys.argv[2]
raw = open(src).read().splitlines()
line={}; xyz={}; elem={}; adj=collections.defaultdict(set); order=[]; h_tmpl=None
for L in raw:
    if L.startswith(("ATOM","HETATM")):
        L=f"{L:<80}"; s=int(L[6:11])
        line[s]=L; elem[s]=L[76:78].strip()
        xyz[s]=[float(L[30:38]),float(L[38:46]),float(L[46:54])]
        order.append(s)
        if elem[s]=="H" and h_tmpl is None: h_tmpl=L
    elif L.startswith("CONECT"):
        f=[L[6+5*k:11+5*k].strip() for k in range(12)]; f=[int(x) for x in f if x]
        for j in f[1:]: adj[f[0]].add(j); adj[j].add(f[0])

nx=next(s for s in order if elem[s]=="N" and any(elem[j]=="O" and len(adj[j])==1 for j in adj[s]))
ox=next(j for j in adj[nx] if elem[j]=="O" and len(adj[j])==1)
ang=lambda a,b,c: math.degrees(math.acos(max(-1,min(1,dot(unit(sub(xyz[a],xyz[b])),unit(sub(xyz[c],xyz[b])))))))
alpha=[j for j in sorted(adj[nx]) if elem[j]=="C"]
print(f"nitroxide N{nx} (formal {line[nx][78:80].strip() or '0'})  O{ox}")
print(f"  in : N-O {nrm(sub(xyz[ox],xyz[nx])):.3f} A, angle sum at N "
      f"{ang(alpha[0],nx,alpha[1])+ang(alpha[0],nx,ox)+ang(alpha[1],nx,ox):.1f} deg")

# 1. promote each alpha C-H to a methyl -------------------------------------
nxt=max(order)+1; new=[]; promoted=[]; protonated=[]
for ca in alpha:
    hs=[j for j in adj[ca] if elem[j]=="H"]
    if len(hs)!=1:
        print(f"  C{ca}: already quaternary ({len(hs)} H), left alone"); continue
    h=hs[0]
    xyz[h]=add(xyz[ca], mul(unit(sub(xyz[h],xyz[ca])), CC))
    elem[h]="C"; promoted.append(h)
    axis=unit(sub(xyz[ca],xyz[h]))
    seed=[1.,0.,0.] if abs(axis[0])<0.9 else [0.,1.,0.]
    u=unit(cross(axis,seed)); v=cross(axis,u)
    for k in range(3):
        phi=k*2*math.pi/3
        d=add(mul(axis,COS_TET), mul(add(mul(u,math.cos(phi)),mul(v,math.sin(phi))),SIN_TET))
        xyz[nxt]=add(xyz[h], mul(unit(d),CH)); elem[nxt]="H"
        adj[h].add(nxt); adj[nxt].add(h); order.append(nxt); new.append(nxt); nxt+=1
    print(f"  C{ca}: H{h} -> methyl carbon, 3 H added")

# 2. protonate the 4-amino linker nitrogen --------------------------------
# The published residue is +1 because this nitrogen is an ammonium. A neutral
# secondary amine here gives a net-0 residue, and no amount of charge
# redistribution substitutes for the missing hydrogen.
amine = [s_ for s_ in order
         if elem[s_] == "N" and s_ != nx
         and not any(elem[j] == "O" for j in adj[s_])]
for n in amine:
    nb = sorted(adj[n])
    hs = [j for j in nb if elem[j] == "H"]
    if len(nb) == 4:
        # Still record it: RDKit needs the formal charge or sanitisation fails
        # on a four-bonded nitrogen, which is what breaks a naive re-run.
        protonated.append(n)
        print(f"  N{n}: already protonated ({len(hs)} H), left alone"); continue
    if len(nb) != 3 or len(hs) != 1:
        print(f"  N{n}: {len(nb)} bonds / {len(hs)} H - not a secondary amine, left alone")
        continue
    # fourth tetrahedral vertex: opposite the sum of the three existing bonds
    d = unit(mul(add(add(unit(sub(xyz[nb[0]], xyz[n])), unit(sub(xyz[nb[1]], xyz[n]))),
                     unit(sub(xyz[nb[2]], xyz[n]))), -1.0))
    xyz[nxt] = add(xyz[n], mul(d, NH)); elem[nxt] = "H"
    adj[n].add(nxt); adj[nxt].add(n); order.append(nxt); new.append(nxt)
    protonated.append(n)
    print(f"  N{n}: H{nxt} added -> ammonium, residue becomes +1")
    nxt += 1

# 3. aminoxyl geometry: planar N, correct N-O ------------------------------
bis=unit(add(unit(sub(xyz[alpha[0]],xyz[nx])), unit(sub(xyz[alpha[1]],xyz[nx]))))
xyz[ox]=add(xyz[nx], mul(bis,-NO))

# 4. MMFF: free the hydrogens and every methyl carbon, fix the rest ---------
idx={}; m=Chem.RWMol()
for s in order:
    a=Chem.Atom(elem[s]); a.SetNoImplicit(True)
    if s in protonated: a.SetFormalCharge(1)      # ammonium N, four bonds
    idx[s]=m.AddAtom(a)
for i in list(adj):
    for j in adj[i]:
        if i<j and i in idx and j in idx: m.AddBond(idx[i],idx[j],Chem.BondType.SINGLE)
cap=m.AddAtom(Chem.Atom("H")); m.AddBond(idx[ox],cap,Chem.BondType.SINGLE)   # MMFF needs a filled valence
for s in order:
    if elem[s]!="C": continue
    os_=[j for j in adj[s] if elem[j]=="O"]
    if len(os_)==2:
        t=[j for j in os_ if not any(elem[k]=="H" for k in adj[j])]
        if t: m.GetBondBetweenAtoms(idx[s],idx[t[0]]).SetBondType(Chem.BondType.DOUBLE)
mol=m.GetMol(); conf=Chem.Conformer(mol.GetNumAtoms())
for s in order: conf.SetAtomPosition(idx[s], Point3D(*xyz[s]))
conf.SetAtomPosition(cap, Point3D(*add(xyz[ox], mul(unit(sub(xyz[ox],xyz[nx])),0.96))))
mol.AddConformer(conf); Chem.SanitizeMol(mol)
print("  smiles:", Chem.MolToSmiles(Chem.RemoveHs(mol)))

free={idx[s] for s in order if elem[s]=="H"} | {cap}
free |= {idx[s] for s in order if elem[s]=="C" and sum(1 for j in adj[s] if elem[j]=="H")==3}
props=AllChem.MMFFGetMoleculeProperties(mol)
ff=AllChem.MMFFGetMoleculeForceField(mol, props)
for a in range(mol.GetNumAtoms()):
    if a not in free: ff.AddFixedPoint(a)
ff.Initialize(); e0=ff.CalcEnergy()
ok=ff.Minimize(maxIts=5000, forceTol=1e-6)==0
print(f"  MMFF: {len(free)} free / {mol.GetNumAtoms()-len(free)} fixed, "
      f"E {e0:.1f} -> {ff.CalcEnergy():.1f} kcal/mol, converged={ok}")
c=mol.GetConformer()
for s in order:
    p=c.GetAtomPosition(idx[s]); xyz[s]=[p.x,p.y,p.z]

# 5. write, splicing only the columns that changed --------------------------
def splice(L,s):
    L=f"{L:<80}"
    if s in protonated:                            # mark the ammonium for CGenFF
        L=L[:78]+"1+"+L[80:]
    L=L[:6]+f"{s:5d}"+L[11:]
    L=L[:12]+f" {elem[s]:<3}"[:4]+L[16:]
    L=L[:30]+f"{xyz[s][0]:8.3f}{xyz[s][1]:8.3f}{xyz[s][2]:8.3f}"+L[54:]
    return (L[:76]+f"{elem[s]:>2}"+L[78:]).rstrip()
out=[L for L in raw if not L.startswith(("ATOM","HETATM","CONECT","END","MASTER"))]
out+=[splice(line[s] if s in line else h_tmpl, s) for s in order]
out+=["CONECT"+f"{s:5d}"+"".join(f"{j:5d}" for j in sorted(adj[s])) for s in order if adj[s]]
out+=["END"]
open(dst,"w").write("\n".join(out)+"\n")

print(f"  out: N-O {nrm(sub(xyz[ox],xyz[nx])):.3f} A, angle sum at N "
      f"{ang(alpha[0],nx,alpha[1])+ang(alpha[0],nx,ox)+ang(alpha[1],nx,ox):.1f} deg")
print(f"  wrote {dst}: {len(order)} atoms ({len(new)} added, {len(promoted)} promoted)")
