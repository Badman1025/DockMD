import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg
from utils import PROJECT_DIR, run, save_state, load_state, get_path

OBABEL = cfg.TOOLS["obabel"]
GMX = cfg.TOOLS["gmx"]
ACPYPE = cfg.TOOLS["acpype"]


def find_latest_docking():
    base = os.path.join(PROJECT_DIR, cfg.RESULT_DIR, cfg.DOCKING_DIR)
    if not os.path.isdir(base):
        sys.exit(f"!! Directory not found: {base}")
    subdirs = [d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d))]
    subdirs.sort(reverse=True)
    return os.path.join(base, subdirs[0])


def extract_molname_from_itp(itp_path):
    with open(itp_path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "[ moleculetype ]" in line:
            for j in range(i+1, min(i+4, len(lines))):
                toks = lines[j].split()
                if toks and not toks[0].startswith(";"):
                    return toks[0]
    return "UNL"


def extract_resname_from_itp(itp_path):
    with open(itp_path) as f:
        in_atoms = False
        for line in f:
            s = line.strip()
            if s.startswith("["):
                in_atoms = s.startswith("[ atoms ]")
                continue
            if in_atoms and s and not s.startswith(";") and not s.startswith("#"):
                toks = s.split()
                if len(toks) > 3:
                    return toks[3]
    return "UNL"


def _safe_int(text, default=None):
    try:
        return int(text)
    except (TypeError, ValueError):
        return default


def _safe_float(text, default=None):
    try:
        return float(text)
    except (TypeError, ValueError):
        return default


def _report_numbering_gaps(pdb_path):
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and len(line) > 25:
                ch = line[21]
                resid = _safe_int(line[22:26])
                if resid is None:
                    continue
                lst = residues.setdefault(ch, [])
                if not lst or lst[-1] != resid:
                    lst.append(resid)
    gaps = []
    for ch in sorted(residues):
        for r1, r2 in zip(residues[ch], residues[ch][1:]):
            if r2 > r1 + 1:
                gaps.append((ch, r1, r2, r2 - r1 - 1))
    if gaps:
        for ch, r1, r2, n in gaps:
            print(f"  [QC] chain {ch}: numbering gap {r1} -> {r2} "
                  f"({n} residue(s) missing in input)")
    else:
        print("  [QC] no numbering gaps found in input chains")


def _check_backbone_closure(pdb_path, cutoff=3.0):
    last_c = {}
    problems = []
    loose = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("TER"):
                last_c.clear()
                continue
            if not line.startswith("ATOM") or len(line) < 55:
                continue
            ch = line[21]
            name = line[12:16].strip()
            resid = _safe_int(line[22:26])
            x = _safe_float(line[30:38])
            y = _safe_float(line[38:46])
            z = _safe_float(line[46:54])
            if resid is None or x is None or y is None or z is None:
                continue
            if name == "C":
                last_c[ch] = (resid, x, y, z)
            elif name == "N" and ch in last_c:
                pr, px, py, pz = last_c[ch]
                d = ((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2) ** 0.5
                if d > cutoff:
                    problems.append((ch, pr, resid, d))
                elif d > 1.8:
                    loose.append((ch, pr, resid, d))
    for ch, r1, r2, d in loose:
        print(f"  [QC] warning: loose joint chain {ch}: C {r1} -> N {r2}: {d:.2f} A")
    if problems:
        print("!! Protein QC failed: unresolved backbone gap(s) after PDBFixer:")
        for ch, r1, r2, d in problems:
            print(f"   chain {ch}: C of residue {r1} -> N of residue {r2}: "
                  f"{d:.2f} A (peptide bond should be ~1.3 A)")
        sys.exit("   PDBFixer could not rebuild the missing residues. "
                 "Provide a structurally complete input PDB and rerun. "
                 "MD will not start on a broken protein.")
    print("  [QC] backbone closure check passed")


def _parse_mol2_atoms(mol2_path):
    atoms = []
    in_atom = False
    try:
        with open(mol2_path) as f:
            for line in f:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atom = True
                    continue
                if line.startswith("@<TRIPOS>"):
                    if in_atom:
                        break
                    continue
                if in_atom:
                    toks = line.split()
                    if len(toks) < 6:
                        continue
                    x = _safe_float(toks[2])
                    y = _safe_float(toks[3])
                    z = _safe_float(toks[4])
                    if x is None or y is None or z is None:
                        continue
                    atoms.append((toks[1], toks[5], x, y, z))
    except OSError:
        return []
    return atoms


def _mol2_element(atom_type):
    el = atom_type.split(".")[0].strip("0123456789").upper()
    return el or "X"


def _parse_mol2_bonds(mol2_path, n_atoms):
    bonds = set()
    in_bond = False
    try:
        with open(mol2_path) as f:
            for line in f:
                if line.startswith("@<TRIPOS>BOND"):
                    in_bond = True
                    continue
                if line.startswith("@<TRIPOS>"):
                    if in_bond:
                        break
                    continue
                if in_bond:
                    toks = line.split()
                    if len(toks) < 4:
                        continue
                    a1 = _safe_int(toks[1])
                    a2 = _safe_int(toks[2])
                    if a1 is None or a2 is None:
                        return None
                    if not (1 <= a1 <= n_atoms and 1 <= a2 <= n_atoms):
                        return None
                    bonds.add((min(a1, a2), max(a1, a2)))
    except OSError:
        return None
    return bonds or None


def _pdb_heavy_coords(pdb_path):
    coords = []
    try:
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                element = line[76:78].strip().upper()
                name = line[12:16].strip()
                if element:
                    is_h = element == "H"
                else:
                    is_h = name.lstrip("0123456789").upper().startswith("H")
                x = _safe_float(line[30:38])
                y = _safe_float(line[38:46])
                z = _safe_float(line[46:54])
                if is_h or x is None or y is None or z is None:
                    continue
                coords.append((round(x, 2), round(y, 2), round(z, 2)))
    except OSError:
        return []
    return coords


def _pdb_element_counts(pdb_path):
    counts = {}
    try:
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                element = line[76:78].strip().upper()
                if not element:
                    element = line[12:16].strip().lstrip("0123456789")[:1].upper()
                if not element or element == "H":
                    continue
                counts[element] = counts.get(element, 0) + 1
    except OSError:
        return {}
    return counts


def _ligand_geometry_qc(mol2_path, ref_pdb):
    try:
        return _ligand_geometry_qc_impl(mol2_path, ref_pdb)
    except Exception as exc:
        return True, f"geometry check unavailable ({exc}); continuing"


def _ligand_geometry_qc_impl(mol2_path, ref_pdb):
    if not os.path.isfile(mol2_path):
        return False, "mol2 file was not generated"
    atoms = _parse_mol2_atoms(mol2_path)
    if not atoms:
        return False, "no ATOM records parsed from mol2"
    bonds = _parse_mol2_bonds(mol2_path, len(atoms))
    nonbonded_cut = 0.9 if bonds else 0.5
    worst = ("", "", 1e9)
    worst_bond = ("", "", 1e9)
    for i in range(len(atoms)):
        _, _, x1, y1, z1 = atoms[i]
        for j in range(i + 1, len(atoms)):
            _, _, x2, y2, z2 = atoms[j]
            d = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5
            if bonds and (i + 1, j + 1) in bonds:
                if d < worst_bond[2]:
                    worst_bond = (atoms[i][0], atoms[j][0], d)
                continue
            if d < worst[2]:
                worst = (atoms[i][0], atoms[j][0], d)
    if worst_bond[2] < 0.8:
        return False, (f"bonded atoms too close: {worst_bond[0]} -- "
                       f"{worst_bond[1]} ({worst_bond[2]:.3f} A, shorter than "
                       f"any real bond)")
    if worst[2] < nonbonded_cut:
        kind = "non-bonded atoms too close" if bonds else "atoms too close"
        return False, (f"{kind}: {worst[0]} -- {worst[1]} ({worst[2]:.3f} A)")

    ref_heavy = sorted(_pdb_heavy_coords(ref_pdb))
    new_heavy = sorted((round(x, 2), round(y, 2), round(z, 2))
                       for _, atype, x, y, z in atoms
                       if _mol2_element(atype) != "H")
    if ref_heavy and len(ref_heavy) != len(new_heavy):
        return False, (f"heavy atom count changed: input {len(ref_heavy)}, "
                       f"mol2 {len(new_heavy)}")
    if ref_heavy:
        for r, n in zip(ref_heavy, new_heavy):
            drift = max(abs(a - b) for a, b in zip(r, n))
            if drift > 0.05:
                return False, (f"heavy atom moved by {drift:.3f} A "
                               f"(docking pose must be preserved)")

    n_h = 0
    new_comp = {}
    for _, atype, _, _, _ in atoms:
        el = _mol2_element(atype)
        if el == "H":
            n_h += 1
        else:
            new_comp[el] = new_comp.get(el, 0) + 1
    print(f"  [QC] ligand: {len(new_heavy)} heavy atoms, {n_h} hydrogens, "
          f"bond block parsed: {'yes' if bonds else 'no'}")
    ref_comp = _pdb_element_counts(ref_pdb)
    if ref_comp and new_comp and ref_comp != new_comp:
        print(f"  [QC] warning: element composition differs from input "
              f"(input {dict(sorted(ref_comp.items()))}, "
              f"mol2 {dict(sorted(new_comp.items()))})")
    if n_h < len(new_heavy) // 2:
        print(f"  [QC] warning: few hydrogens for this ligand "
              f"({n_h} H on {len(new_heavy)} heavy atoms)")
    return True, (f"{len(atoms)} atoms, min non-bonded distance "
                  f"{worst[2]:.2f} A, heavy atoms unchanged")


def _report_interface_contact(md_dir, n_prot, n_lig):
    try:
        path = os.path.join(md_dir, "complex_initial.gro")
        if not os.path.isfile(path):
            return
        with open(path) as f:
            lines = f.readlines()
        declared = _safe_int(lines[1].strip())
        if declared is None or declared != n_prot + n_lig:
            return
        coords = []
        for line in lines[2:2 + declared]:
            if len(line) < 44:
                return
            name = line[10:15].strip()
            x = _safe_float(line[20:28].strip())
            y = _safe_float(line[28:36].strip())
            z = _safe_float(line[36:44].strip())
            if x is None or y is None or z is None:
                return
            coords.append((name, x, y, z))
        prot = [(x, y, z) for name, x, y, z in coords[:n_prot]
                if not name.upper().startswith("H")]
        lig = [(x, y, z) for name, x, y, z in coords[n_prot:]
               if not name.upper().startswith("H")]
        if not prot or not lig:
            return
        if len(prot) * len(lig) > 40000000:
            return
        best = 1e9
        for lx, ly, lz in lig:
            for px, py, pz in prot:
                d = ((lx - px) ** 2 + (ly - py) ** 2 + (lz - pz) ** 2) ** 0.5
                if d < best:
                    best = d
        if best < 1.5:
            print(f"  [QC] warning: ligand-protein closest heavy-atom contact "
                  f"{best:.2f} A (severe overlap; inspect final_check.pdb)")
        elif best < 2.0:
            print(f"  [QC] note: ligand-protein closest heavy-atom contact "
                  f"{best:.2f} A (tight, minimization should relieve it)")
        else:
            print(f"  [QC] ligand-protein closest heavy-atom contact "
                  f"{best:.2f} A")
    except Exception as exc:
        print(f"  [QC] interface contact check skipped ({exc})")


def strip_atomtypes(itp_in, itp_out):
    with open(itp_in) as f:
        lines = f.readlines()
    atomtypes, clean = [], []
    inside = False
    for line in lines:
        s = line.strip()
        if s.startswith("[ atomtypes ]"):
            inside = True
            atomtypes.append(line)
            continue
        if inside and s.startswith("[") and not s.startswith("[ atomtypes ]"):
            inside = False
        if inside:
            atomtypes.append(line)
        else:
            clean.append(line)
    with open(itp_out, "w") as f:
        f.writelines(clean)
    return atomtypes


def _try_obabel_hydrogens(md_dir):
    try:
        run([OBABEL, "ligand.pdb", "-O", "lig_noh.mol2", "-d"], cwd=md_dir)
        run([OBABEL, "lig_noh.mol2", "-O", "lig.mol2", "-h"], cwd=md_dir)
    except SystemExit:
        return False, "obabel conversion failed"
    return _ligand_geometry_qc(os.path.join(md_dir, "lig.mol2"),
                               os.path.join(md_dir, "ligand.pdb"))


def _try_pymol_hydrogens(md_dir, lig_mol2):
    script = os.path.join(md_dir, "_lig_h.py")
    try:
        with open(script, "w") as f:
            f.write("from pymol import cmd\n")
            f.write(f'cmd.load(r"{os.path.join(md_dir, "ligand.pdb")}", "lig")\n')
            f.write('cmd.remove("hydrogens")\n')
            f.write('cmd.h_add("lig")\n')
            f.write(f'cmd.save(r"{lig_mol2}", "lig")\n')
        run([cfg.TOOLS["pymol"], "-cq", script], cwd=md_dir)
    except SystemExit:
        return False, "pymol hydrogen placement failed"
    finally:
        try:
            os.remove(script)
        except OSError:
            pass
    return _ligand_geometry_qc(lig_mol2, os.path.join(md_dir, "ligand.pdb"))


def main():
    chains = "".join(cfg.PROTEIN_CHAINS.split()).upper()
    if not chains:
        sys.exit("!! PROTEIN_CHAINS in config.py is empty (e.g. 'A' or 'AB')")
    print(f"[0] Target chains: {'+'.join(chains)} ({len(chains)} chain(s))")

    docking_dir = get_path(load_state(), "docking_dir",
                           fallback_finder=find_latest_docking,
                           what="Docking result directory")
    print(f"[1] Using docking result directory: {docking_dir}")

    lig_pdb = os.path.join(docking_dir, f"{cfg.LIGAND_NAME}_docked.pdb")
    if not os.path.isfile(lig_pdb):
        sys.exit(f"!! Missing file: {lig_pdb}")

    md_dir = os.path.join(PROJECT_DIR, cfg.RESULT_DIR, cfg.MD_SETUP_DIR,
                          f"{cfg.PROTEIN_NAME}_{chains}_{cfg.LIGAND_NAME}")
    os.makedirs(md_dir, exist_ok=True)
    print(f"[2] MD working directory: {md_dir}")

    print(f"\n[3a] Extracting pure chains {'+'.join(chains)} "
          f"(TER between chains, SEQRES preserved)...")
    raw_pdb = os.path.join(PROJECT_DIR, cfg.RAW_DATA_DIR, cfg.PROTEIN_PDB)
    if not os.path.isfile(raw_pdb):
        sys.exit(f"!! Raw protein file not found: {raw_pdb}")
    pure_pdb = os.path.join(md_dir, f"receptor_{chains}.pdb")
    with open(raw_pdb) as f, open(pure_pdb, "w") as g:
        prev_chain = None
        for line in f:
            record = line[:6].strip()
            if record == "SEQRES":
                if len(line) > 11 and line[11] in chains:
                    g.write(line)
            elif record == "ATOM" and len(line) > 21 and line[21] in chains:
                if prev_chain is not None and line[21] != prev_chain:
                    g.write("TER\n")
                g.write(line)
                prev_chain = line[21]
        if prev_chain is not None:
            g.write("TER\n")
        g.write("END\n")
    if prev_chain is None:
        sys.exit(f"!! No ATOM records found for chains {'+'.join(chains)} in {raw_pdb}")

    _report_numbering_gaps(pure_pdb)

    print("\n[3b] Fixing missing atoms with PDBFixer ...")
    fix_script = f"""
from pdbfixer import PDBFixer
from openmm.app import PDBFile
fixer = PDBFixer(filename=r"{pure_pdb}")
fixer.findMissingResidues()
chain_lengths = {{}}
for chain in fixer.topology.chains():
    chain_lengths[chain.index] = sum(1 for _ in chain.residues())
fixer.missingResidues = {{
    key: val for key, val in fixer.missingResidues.items()
    if 0 < key[1] < chain_lengths[key[0]]
}}
fixer.findMissingAtoms()
fixer.addMissingAtoms()
with open(r"{os.path.join(md_dir, 'receptor_fixed.pdb')}", "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)
"""
    tmp_fix = os.path.join(md_dir, "_fix.py")
    with open(tmp_fix, "w") as f:
        f.write(fix_script)
    run([sys.executable, tmp_fix], cwd=md_dir)
    os.remove(tmp_fix)
    _check_backbone_closure(os.path.join(md_dir, "receptor_fixed.pdb"))

    print("\n[3c] gmx pdb2gmx: building protein topology ...")
    run([GMX, "pdb2gmx", "-f", "receptor_fixed.pdb",
         "-o", "protein.gro", "-p", "topol.top",
         "-i", "posre.itp", "-water", "tip3p",
         "-ff", "amber99sb-ildn", "-ignh"], cwd=md_dir)

    print("\n[4] Ligand parameterization: hydrogen placement + acpype ...")
    shutil.copy(lig_pdb, os.path.join(md_dir, "ligand.pdb"))
    lig_mol2 = os.path.join(md_dir, "lig.mol2")

    print("  [4.1] Strategy 1: obabel (delete hydrogens, then re-add) ...")
    ok, msg = _try_obabel_hydrogens(md_dir)
    if ok:
        print(f"  [QC] ligand geometry passed (obabel): {msg}")
    else:
        print(f"  [QC] strategy 1 rejected: {msg}")
        print("  [4.2] Strategy 2: PyMOL h_add (independent hydrogen "
              "placement algorithm) ...")
        ok, msg = _try_pymol_hydrogens(md_dir, lig_mol2)
        if ok:
            print(f"  [QC] ligand geometry passed (pymol): {msg}")
        else:
            print(f"  [QC] strategy 2 rejected: {msg}")
            sys.exit(
                f"!! Ligand QC failed for {cfg.LIGAND_NAME}: obabel and PyMOL "
                f"both produced an invalid mol2 geometry.\n"
                f"   The docking pose itself is geometrically clean, so this "
                f"is a hydrogen-placement limitation of both converters.\n"
                f"   Options: (1) provide a valid mol2 for this ligand and "
                f"rerun prep with it in the docking folder;\n"
                f"            (2) check the protonation state of the input "
                f"sdf in 1.Raw_data.\n"
                f"   MD will not start on a broken ligand.")

    run([ACPYPE, "-i", "lig.mol2", "-b", cfg.LIGAND_NAME,
         "-n", "0", "-a", "gaff2"], cwd=md_dir)

    acpype_dir = os.path.join(md_dir, f"{cfg.LIGAND_NAME}.acpype")
    itp_full = os.path.join(acpype_dir, f"{cfg.LIGAND_NAME}_GMX.itp")
    itp_clean = os.path.join(acpype_dir, f"{cfg.LIGAND_NAME}_clean.itp")
    gro_lig = os.path.join(acpype_dir, f"{cfg.LIGAND_NAME}_GMX.gro")

    if not os.path.isfile(itp_full):
        sys.exit(f"!! acpype output missing: {itp_full}")

    mol_name = extract_molname_from_itp(itp_full)
    print(f"  Ligand molecule name: {mol_name}")

    lig_res = extract_resname_from_itp(itp_full)
    print(f"  Ligand residue name: {lig_res}")

    atomtypes_block = strip_atomtypes(itp_full, itp_clean)
    print(f"  Stripped {len(atomtypes_block)} atomtypes lines from itp")

    print("\n[5] Concatenating protein.gro + ligand gro ...")
    protein_gro = os.path.join(md_dir, "protein.gro")
    if not os.path.isfile(protein_gro) or not os.path.isfile(gro_lig):
        sys.exit("!! Missing protein.gro or ligand gro from acpype")
    with open(protein_gro) as f:
        p_lines = f.readlines()
    with open(gro_lig) as f:
        l_lines = f.readlines()

    n_p = _safe_int(p_lines[1].strip()) if len(p_lines) > 1 else None
    n_l = _safe_int(l_lines[1].strip()) if len(l_lines) > 1 else None
    if not n_p or not n_l:
        sys.exit("!! Could not read atom counts from protein.gro / ligand gro")
    if len(p_lines) < 3 + n_p or len(l_lines) < 3 + n_l:
        sys.exit("!! Truncated gro file: atom count does not match content")

    with open(os.path.join(md_dir, "complex_initial.gro"), "w") as f:
        f.write("Protein-Ligand Complex\n")
        f.write(f"{n_p + n_l}\n")
        f.writelines(p_lines[2:2 + n_p])
        f.writelines(l_lines[2:2 + n_l])
        f.write(p_lines[2 + n_p])

    _report_interface_contact(md_dir, n_p, n_l)

    print("\n[6] Updating topol.top ...")
    with open(os.path.join(md_dir, "topol.top")) as f:
        top = f.readlines()

    protein_mol_lines = []
    in_mols = False
    for line in top:
        s = line.strip()
        if s.startswith("[ molecules ]"):
            in_mols = True
            continue
        if in_mols:
            if s.startswith("["):
                break
            if s and not s.startswith(";") and not s.startswith("#"):
                protein_mol_lines.append(line)
    if not protein_mol_lines:
        sys.exit("!! Could not find protein entries in topol.top [ molecules ]")
    print(f"  Protein molecule entries from pdb2gmx: "
          f"{[l.split()[0] for l in protein_mol_lines]}")

    new_top = []
    inserted_atomtypes = False
    inserted_include = False
    for line in top:
        if not inserted_atomtypes and line.strip().startswith('#include "amber99sb-ildn.ff/forcefield.itp"'):
            new_top.append(line)
            new_top.append("\n; === Ligand atomtypes ===\n")
            new_top.extend(atomtypes_block)
            new_top.append("\n")
            inserted_atomtypes = True
            continue
        if not inserted_include and line.strip().startswith("[ system ]"):
            new_top.append(f'; Include ligand topology\n#include "{cfg.LIGAND_NAME}.acpype/{cfg.LIGAND_NAME}_clean.itp"\n\n')
            inserted_include = True
        if line.strip().startswith("[ molecules ]"):
            new_top.append("[ molecules ]\n")
            new_top.append("; Compound        #mols\n")
            new_top.extend(protein_mol_lines)
            new_top.append(f"{mol_name:<20} 1\n")
            break
        new_top.append(line)

    with open(os.path.join(md_dir, "topol.top"), "w") as f:
        f.writelines(new_top)

    print("\n[7] editconf: building box ...")
    run([GMX, "editconf", "-f", "complex_initial.gro",
         "-o", "complex_box.gro", "-bt", "cubic", "-d", "1.0", "-c"], cwd=md_dir)

    print("\n[8] solvate: adding water ...")
    run([GMX, "solvate", "-cp", "complex_box.gro",
         "-cs", "spc216.gro", "-o", "complex_solv.gro",
         "-p", "topol.top"], cwd=md_dir)

    print("\n[9] Generating final_check.pdb for PyMOL inspection ...")
    run([GMX, "editconf", "-f", "complex_solv.gro",
         "-o", "final_check.pdb"], cwd=md_dir)

    save_state(md_dir=md_dir, ligand_name=cfg.LIGAND_NAME,
               lig_resname=lig_res, protein_name=cfg.PROTEIN_NAME,
               protein_chains=chains)
    print(f"\n[OK] All done! MD directory: {md_dir}")


if __name__ == "__main__":
    main()
