import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg
from utils import PROJECT_DIR

WATER_RESNAMES = frozenset({"HOH", "WAT", "DOD", "H2O", "SOL", "TIP", "TIP3"})
TWO_LETTER_ELEMENTS = frozenset({
    "LI", "NA", "MG", "AL", "SI", "CL", "CA", "MN", "FE", "CO", "NI",
    "CU", "ZN", "SE", "BR", "CD", "HG", "CS", "BA", "PT", "AU", "AG"})
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}
CONTACT_CUTOFF = 5.0
INTERFACE_CUTOFF = 4.5
INTERFACE_SEARCH = 8.0
DISULFIDE_CUTOFF = 2.5
BOX_MARGIN = 6.0
MIN_BOX = 20.0


def _element_from_pdb_name(name, is_het):
    letters = ""
    for ch in name:
        if ch.isalpha():
            letters += ch
        elif letters:
            break
    if not letters:
        return "X"
    if is_het and len(letters) >= 2 and letters[:2].upper() in TWO_LETTER_ELEMENTS:
        return letters[:2].upper()
    return letters[0].upper()


def _element_from_mol2_type(atom_type):
    head = atom_type.split(".")[0]
    if not head:
        return "X"
    if len(head) >= 2 and head[1].islower():
        return head[:2].upper()
    return head[:1].upper()


def read_pdb(path):
    atoms = []
    residues = []
    seen_res = set()
    skipped = 0
    seen_model = False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "MODEL":
                if seen_model:
                    break
                seen_model = True
                continue
            if rec == "ENDMDL" and seen_model:
                break
            if rec not in ("ATOM", "HETATM"):
                continue
            if len(line) > 80 and line[16] not in (" ", "A"):
                skipped += 1
                continue
            is_het = rec == "HETATM"
            name = line[12:16].strip()
            element = line[76:78].strip().upper() if len(line) > 78 else ""
            if not element:
                element = _element_from_pdb_name(name, is_het)
            chain = line[21] if len(line) > 21 else ""
            chain = chain if chain.strip() else ""
            icode = line[26].strip() if len(line) > 26 else ""
            resseq = line[22:26].strip()
            atoms.append({
                "het": is_het,
                "name": name,
                "resname": line[17:20].strip().upper(),
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "element": element,
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            })
            if not is_het:
                key = (chain, resseq, icode)
                if key not in seen_res:
                    seen_res.add(key)
                    residues.append(key)
    return atoms, skipped, residues


def read_sdf(path):
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    if len(lines) < 4:
        sys.exit(f"!! Ligand file is too short: {path}")
    atoms = []
    if "V3000" in lines[3]:
        in_atom = False
        for line in lines:
            if "BEGIN ATOM" in line:
                in_atom = True
                continue
            if "END ATOM" in line:
                in_atom = False
                continue
            if not in_atom or not line.startswith("M  V30 "):
                continue
            parts = line[7:].split()
            if len(parts) < 6 or not parts[0].isdigit():
                continue
            atoms.append({
                "element": parts[1].upper(),
                "xyz": (float(parts[2]), float(parts[3]), float(parts[4])),
                "name": parts[1] + parts[0],
            })
        return atoms
    try:
        n_atoms = int(lines[3][:3])
    except ValueError:
        sys.exit(f"!! Unrecognized SDF counts line in {path}: {lines[3]!r}")
    for i in range(4, min(4 + n_atoms, len(lines))):
        parts = lines[i].split()
        if len(parts) < 4:
            continue
        atoms.append({
            "element": parts[3].upper(),
            "xyz": (float(parts[0]), float(parts[1]), float(parts[2])),
            "name": parts[3] + str(len(atoms) + 1),
        })
    return atoms


def read_mol2(path):
    atoms = []
    in_atom = False
    for line in open(path, encoding="utf-8", errors="replace"):
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            in_atom = s == "@<TRIPOS>ATOM"
            continue
        if not in_atom or not s:
            continue
        parts = s.split()
        if len(parts) < 6:
            continue
        atoms.append({
            "element": _element_from_mol2_type(parts[5]),
            "xyz": (float(parts[2]), float(parts[3]), float(parts[4])),
            "name": parts[1],
        })
    return atoms


def read_ligand(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".sdf", ".mol"):
        atoms = read_sdf(path)
    elif ext == ".mol2":
        atoms = read_mol2(path)
    elif ext in (".pdb", ".ent"):
        atoms, _skipped, _res = read_pdb(path)
    else:
        sys.exit(f"!! Unsupported ligand format: {ext} (use .sdf, .mol, .mol2 or .pdb)")
    if not atoms:
        sys.exit(f"!! No atoms parsed from ligand file: {path}")
    return atoms


class Grid:
    def __init__(self, points, cell):
        self.cell = cell
        self.points = points
        self.buckets = {}
        for idx, p in enumerate(points):
            self.buckets.setdefault(self._key(p), []).append(idx)

    def _key(self, p):
        c = self.cell
        return (int(math.floor(p[0] / c)), int(math.floor(p[1] / c)), int(math.floor(p[2] / c)))

    def near(self, p, radius):
        c = self.cell
        span = int(math.ceil(radius / c))
        kx, ky, kz = self._key(p)
        r2 = radius * radius
        out = []
        for dx in range(-span, span + 1):
            for dy in range(-span, span + 1):
                for dz in range(-span, span + 1):
                    for idx in self.buckets.get((kx + dx, ky + dy, kz + dz), ()):
                        q = self.points[idx]
                        d2 = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
                        if d2 <= r2:
                            out.append((idx, math.sqrt(d2)))
        return out


def centroid(points):
    n = len(points)
    return tuple(sum(p[i] for p in points) / n for i in range(3))


def extents(points):
    return tuple(max(p[i] for p in points) - min(p[i] for p in points) for i in range(3))


def element_counts(atoms):
    counts = {}
    for a in atoms:
        counts[a["element"]] = counts.get(a["element"], 0) + 1
    return counts


def chain_inventory(atoms, residues):
    order = []
    stats = {}
    for a in atoms:
        if a["het"]:
            continue
        rec = stats.get(a["chain"])
        if rec is None:
            rec = {"atoms": 0, "heavy": 0, "h": 0, "res": []}
            stats[a["chain"]] = rec
            order.append(a["chain"])
        rec["atoms"] += 1
        if a["element"] == "H":
            rec["h"] += 1
        else:
            rec["heavy"] += 1
    for chain, resseq, icode in residues:
        if chain in stats:
            stats[chain]["res"].append(resseq + icode)
    for chain in order:
        rec = stats[chain]
        rec["n_res"] = len(rec["res"])
        rec["range"] = f"{rec['res'][0]}-{rec['res'][-1]}" if rec["res"] else "n/a"
    return order, stats


def residue_sequence(atoms, chain):
    seq = []
    seen = set()
    for a in atoms:
        if a["het"] or a["chain"] != chain:
            continue
        key = (a["resseq"], a["icode"])
        if key in seen:
            continue
        seen.add(key)
        seq.append(AA3_TO_1.get(a["resname"], "X"))
    return "".join(seq)


def write_inventory_csv(out_path, order, inventory, atoms):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["chain", "residues", "heavy_atoms", "hydrogens",
                    "first_residue", "last_residue", "sequence"])
        for chain in order:
            rec = inventory[chain]
            seq = residue_sequence(atoms, chain)
            first = rec["res"][0] if rec["res"] else ""
            last = rec["res"][-1] if rec["res"] else ""
            w.writerow([chain or "-", rec["n_res"], rec["heavy"], rec["h"],
                        first, last, seq])
    return out_path


def het_groups(atoms):
    groups = {}
    for a in atoms:
        if not a["het"]:
            continue
        groups.setdefault((a["chain"], a["resname"], a["resseq"]), []).append(a)
    return groups


def heavy_by_chain(atoms):
    out = {}
    for a in atoms:
        if a["het"] or a["element"] == "H":
            continue
        out.setdefault(a["chain"], []).append(a)
    return out


def pocket_stats(lig_xyz, chain_atoms):
    points = [a["xyz"] for a in chain_atoms]
    grid = Grid(points, CONTACT_CUTOFF)
    hit_atoms = set()
    residues = {}
    for p in lig_xyz:
        for idx, d in grid.near(p, CONTACT_CUTOFF):
            hit_atoms.add(idx)
            key = (chain_atoms[idx]["resname"], chain_atoms[idx]["resseq"])
            if key not in residues or d < residues[key]:
                residues[key] = d
    min_dist = min(residues.values()) if residues else float("inf")
    if min_dist == float("inf"):
        for p in lig_xyz:
            for q in points:
                d = math.dist(p, q)
                if d < min_dist:
                    min_dist = d
            if min_dist <= 1.5:
                break
    return {
        "n_atoms": len(hit_atoms),
        "n_res": len(residues),
        "min_dist": min_dist,
        "residues": sorted(residues.items(), key=lambda kv: kv[1]),
    }


def disulfides(atoms):
    sg = [a for a in atoms if a["element"] == "S" and a["name"].upper() == "SG"]
    bonds = []
    used = set()
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            if i in used or j in used:
                continue
            d = math.dist(sg[i]["xyz"], sg[j]["xyz"])
            if d <= DISULFIDE_CUTOFF:
                bonds.append((sg[i], sg[j], d))
                used.add(i)
                used.add(j)
    free = [a for i, a in enumerate(sg) if i not in used]
    return bonds, free


def chain_interfaces(chain_atoms, order):
    pairs = []
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            ca = chain_atoms.get(order[i], [])
            cb = chain_atoms.get(order[j], [])
            if not ca or not cb:
                continue
            grid = Grid([a["xyz"] for a in cb], INTERFACE_SEARCH)
            best = float("inf")
            n_close = 0
            for a in ca:
                for _idx, d in grid.near(a["xyz"], INTERFACE_SEARCH):
                    if d < best:
                        best = d
                    if d <= INTERFACE_CUTOFF:
                        n_close += 1
            if best <= INTERFACE_SEARCH:
                pairs.append((order[i], order[j], best, n_close))
    pairs.sort(key=lambda t: t[2])
    return pairs


def format_residues(items, limit=6):
    shown = [f"{name}{num}" for (name, num), _d in items[:limit]]
    if len(items) > limit:
        shown.append(f"+{len(items) - limit} more")
    return ", ".join(shown) if shown else "-"


def resolve_protein(pdb_path):
    raw_dir = os.path.join(PROJECT_DIR, cfg.RAW_DATA_DIR)
    path = pdb_path or os.path.join(raw_dir, cfg.PROTEIN_PDB)
    if not os.path.isfile(path):
        listing = ""
        if os.path.isdir(raw_dir):
            names = sorted(n for n in os.listdir(raw_dir) if n.lower().endswith((".pdb", ".ent")))
            listing = "\n   Available in {}: {}".format(cfg.RAW_DATA_DIR, ", ".join(names) or "none")
        sys.exit(f"!! Protein file not found: {path}{listing}")
    return path


def resolve_ligand(ligand_path, protein_path):
    if ligand_path:
        if not os.path.isfile(ligand_path):
            sys.exit(f"!! Ligand file not found: {ligand_path}")
        return ligand_path
    raw_dir = os.path.join(PROJECT_DIR, cfg.RAW_DATA_DIR)
    configured = os.path.join(raw_dir, cfg.PROTEIN_PDB)
    if os.path.abspath(protein_path) != os.path.abspath(configured):
        return _sibling_ligand(protein_path)
    preferred = os.path.join(raw_dir, f"{cfg.LIGAND_NAME}.sdf")
    if os.path.isfile(preferred):
        return preferred
    if os.path.isdir(raw_dir):
        sdfs = sorted(n for n in os.listdir(raw_dir) if n.lower().endswith(".sdf"))
        if sdfs:
            return os.path.join(raw_dir, sdfs[0])
    return None


def _sibling_ligand(protein_path):
    folder = os.path.dirname(os.path.abspath(protein_path))
    stem = os.path.splitext(os.path.basename(protein_path))[0]
    exact = [os.path.join(folder, stem + ext) for ext in (".sdf", ".mol2", ".mol")]
    for path in exact:
        if os.path.isfile(path):
            return path
    prefix = stem.split("-")[0].split("_")[0]
    if len(prefix) < 3:
        return None
    if os.path.isdir(folder):
        for name in sorted(os.listdir(folder)):
            low = name.lower()
            if not low.endswith((".sdf", ".mol2", ".mol")):
                continue
            if name.startswith(prefix + "-") or name.startswith(prefix + "_"):
                return os.path.join(folder, name)
    return None


def het_extras(groups, lig_xyz):
    lig_center = centroid(lig_xyz)
    extras = set()
    for key, group in groups.items():
        if key[1] in WATER_RESNAMES:
            continue
        center = centroid([a["xyz"] for a in group])
        if math.dist(center, lig_center) < 1.5:
            continue
        extras.add(f"{key[1]}({key[0] or '-'}{key[2]})")
    return sorted(extras)


def main(pdb_path=None, ligand_path=None, ligand_chain=None):
    protein_path = resolve_protein(pdb_path)
    ligand_chain = ligand_chain.strip().upper() if ligand_chain else None
    configured_path = os.path.join(PROJECT_DIR, cfg.RAW_DATA_DIR, cfg.PROTEIN_PDB)
    is_configured = os.path.abspath(protein_path) == os.path.abspath(configured_path)
    atoms, skipped_altloc, residues = read_pdb(protein_path)
    if not atoms:
        sys.exit(f"!! No ATOM/HETATM records found in {protein_path}")

    order, inventory = chain_inventory(atoms, residues)
    groups = het_groups(atoms)
    chain_atoms = heavy_by_chain(atoms)
    warnings = []

    print(f"\n[0] Structure: {os.path.basename(protein_path)}")
    n_nonwater = sum(1 for k in groups if k[1] not in WATER_RESNAMES)
    print(f"    atoms {len(atoms)} | polymer chains {len(order)} "
          f"({' '.join(c or '-' for c in order)}) | HETATM groups {len(groups)} "
          f"({n_nonwater} non-water)")
    h_total = sum(1 for a in atoms if a["element"] == "H")
    if h_total:
        print(f"    hydrogens present: {h_total} (pdb2gmx -ignh rebuilds them)")
    if skipped_altloc:
        print(f"    alternate locations skipped: {skipped_altloc}")

    print("\n[1] Chain inventory")
    print("    chain  residues  heavy  range      sequence")
    for chain in order:
        rec = inventory[chain]
        seq = residue_sequence(atoms, chain)
        if len(seq) > 46:
            seq = seq[:43] + "..."
        tag = "  <- ligand chain" if chain == ligand_chain else ""
        print(f"    {chain or '-':<6} {rec['n_res']:>8} {rec['heavy']:>6}  "
              f"{rec['range']:<9}  {seq}{tag}")

    stem = os.path.splitext(os.path.basename(protein_path))[0]
    csv_path = os.path.join(PROJECT_DIR, cfg.RESULT_DIR,
                            f"{stem}_chain_inventory.csv")
    write_inventory_csv(csv_path, order, inventory, atoms)
    print(f"    saved: {os.path.relpath(csv_path, PROJECT_DIR)}")

    ligand_atoms = None
    ligand_label = ""
    chosen_group = None
    stats = {}
    if ligand_chain:
        target = ligand_chain.strip().upper()
        if target not in order:
            sys.exit(f"!! Chain '{ligand_chain}' not present. "
                     f"Chain IDs in {os.path.basename(protein_path)}: "
                     + "+".join(c or "-" for c in order))
        ligand_atoms = [a for a in atoms if a["chain"] == target]
        ligand_label = f"polymer chain {target} of the input structure"
        analysis_order = [c for c in order if c != target]
    else:
        analysis_order = list(order)
        resolved = resolve_ligand(ligand_path, protein_path)
        if resolved:
            ligand_atoms = read_ligand(resolved)
            ligand_label = os.path.basename(resolved)
        else:
            candidates = [(k, v) for k, v in groups.items() if k[1] not in WATER_RESNAMES]
            if candidates:
                chosen_group, ligand_atoms = max(candidates, key=lambda kv: len(kv[1]))
                ligand_label = (f"internal HETATM {chosen_group[1]} "
                                f"{chosen_group[0] or '-'}{chosen_group[2]}")

    if ligand_atoms is None:
        print("\n[2] Ligand: none found (no --ligand, no *.sdf, no non-water HETATM)")
        warnings.append("no ligand available: FIXED_CENTER cannot be derived, docking would "
                        "fall back to the protein centroid (blind docking)")
    else:
        heavy = [a for a in ligand_atoms if a["element"] != "H"]
        if not heavy:
            sys.exit("!! Ligand contains no heavy atoms")
        lig_xyz = [a["xyz"] for a in heavy]
        cx = centroid(lig_xyz)
        ext = extents(lig_xyz)
        formula = " ".join(f"{el}{n}" for el, n in sorted(element_counts(heavy).items()))
        print(f"\n[2] Ligand: {ligand_label}")
        print(f"    heavy atoms {len(heavy)} | composition {formula}")
        print(f"    centroid ({cx[0]:.3f}, {cx[1]:.3f}, {cx[2]:.3f})")
        print(f"    extent {ext[0]:.1f} x {ext[1]:.1f} x {ext[2]:.1f} A")
        if len(heavy) < 4:
            warnings.append("ligand has fewer than 4 heavy atoms: pocket assignment may be unreliable")
        others = het_extras(groups, lig_xyz)
        if others:
            warnings.append("additional non-water HETATM groups present: " + ", ".join(others))

        print(f"\n[3] Pocket assignment (heavy-atom cutoff {CONTACT_CUTOFF:.1f} A)")
        print("    chain  min_dist  atoms<cut  residues")
        for chain in analysis_order:
            ca = chain_atoms.get(chain, [])
            if not ca:
                continue
            st = pocket_stats(lig_xyz, ca)
            stats[chain] = st
            md = f"{st['min_dist']:.2f}" if st["min_dist"] != float("inf") else "n/a"
            print(f"    {chain or '-':<6} {md:>8} {st['n_atoms']:>10}  {format_residues(st['residues'])}")

    print(f"\n[4] Chain interfaces (heavy-atom cutoff {INTERFACE_CUTOFF:.1f} A)")
    pairs = []
    if len(order) < 2:
        print("    single polymer chain, no interface")
    else:
        pairs = chain_interfaces(chain_atoms, order)
        if not pairs:
            print(f"    no chain pair closer than {INTERFACE_SEARCH:.1f} A")
        for a, b, d, n_close in pairs[:10]:
            print(f"    {a or '-'}-{b or '-'}: min {d:.2f} A, {n_close} atom pairs within cutoff")

    if ligand_atoms is None and pairs:
        together = {p[0] for p in pairs} | {p[1] for p in pairs}
        short = [c for c in order if inventory[c]["n_res"] <= 30 and c in together]
        if short:
            warnings.append("short chain(s) " + "+".join(c or "-" for c in short)
                            + " contact the receptor and may be the ligand: "
                            + f"rerun with --ligand-chain {short[0]}")

    bonds, free_sg = disulfides([a for a in atoms if not a["het"]])
    print(f"\n[5] Disulfides (SG-SG cutoff {DISULFIDE_CUTOFF:.1f} A)")
    if not bonds:
        print("    none detected")
    for a, b, d in bonds:
        kind = "intra" if a["chain"] == b["chain"] else "INTER"
        print(f"    {kind:<5} {a['chain'] or '-'}{a['resseq']} - "
              f"{b['chain'] or '-'}{b['resseq']}  ({d:.2f} A)")
    if free_sg:
        print(f"    cysteines without SG partner: {len(free_sg)}")

    print("\n[6] Suggested config")
    suggested = ""
    if ligand_atoms is not None and stats:
        pocket = [c for c, st in stats.items() if st["n_atoms"] >= 3 or st["min_dist"] <= 3.5]
        if not pocket:
            closest = min(stats.items(), key=lambda kv: kv[1]["min_dist"])
            pocket = [closest[0]]
            warnings.append(f"ligand is farther than {CONTACT_CUTOFF:.0f} A from every chain "
                            f"(closest {closest[0] or '-'} at {closest[1]['min_dist']:.2f} A): "
                            f"check that the ligand file matches the structure")
        keep = list(pocket)
        for a, b, _d in bonds:
            for held, dropped in ((a, b), (b, a)):
                if held["chain"] in keep and dropped["chain"] not in keep:
                    warnings.append(f"chain {dropped['chain'] or '-'} is disulfide-linked to "
                                    f"{held['chain'] or '-'}: removing it leaves a free cysteine")
                    keep.append(dropped["chain"])
        suggested = "".join(c for c in analysis_order if c in keep)
        excluded = [c for c in analysis_order if c not in keep]
        if excluded:
            warnings.append("chain(s) " + "+".join(c or "-" for c in excluded)
                            + " do not contact the ligand and are omitted")
            biggest = max((inventory[c]["n_res"] for c in keep if c in inventory), default=0)
            full = [c for c in excluded
                    if inventory.get(c, {}).get("n_res", 0) >= 0.5 * biggest and biggest >= 40]
            if full:
                warnings.append("omitted chain(s) " + "+".join(c or "-" for c in full)
                                + " are full-length copies, not peptides: "
                                + "they may host an independent binding site")
        if len(pocket) > 1:
            warnings.append("binding site spans " + "+".join(c or "-" for c in pocket)
                            + ": keep all of them, otherwise the pocket is incomplete")

    if suggested:
        print(f'    PROTEIN_CHAINS  = "{suggested}"')
        if is_configured and suggested != cfg.PROTEIN_CHAINS:
            print(f'    note: config.py currently has "{cfg.PROTEIN_CHAINS}"')
    else:
        print("    no suggestion available: select chains manually and set FIXED_CENTER "
              "from a known binding site")

    if ligand_atoms is not None:
        heavy = [a for a in ligand_atoms if a["element"] != "H"]
        cx = centroid([a["xyz"] for a in heavy])
        ext = extents([a["xyz"] for a in heavy])
        box = max(MIN_BOX, math.ceil((max(ext) + 2 * BOX_MARGIN) / 2.0) * 2.0)
        print(f"    FIXED_CENTER    = ({cx[0]:.3f}, {cx[1]:.3f}, {cx[2]:.3f})")
        print(f"    BOX_SIZE        = {int(box)}")
        if is_configured:
            if cfg.FIXED_CENTER is None:
                print("    note: config.py leaves FIXED_CENTER unset (protein-centroid fallback)")
            elif math.dist(cx, cfg.FIXED_CENTER) > 2.0:
                print(f"    note: config.py center is {math.dist(cx, cfg.FIXED_CENTER):.2f} A away "
                      f"from the ligand centroid")
        if box > 40:
            warnings.append(f"ligand is large (max extent {max(ext):.1f} A): "
                            f"a {int(box)} A box will slow the search")

    if any(not c for c in order):
        warnings.append("structure contains blank chain IDs: assign chain identifiers "
                        "before using PROTEIN_CHAINS")

    nonstandard = sorted({a["resname"] for a in atoms
                          if not a["het"] and a["resname"] not in AA3_TO_1})
    if nonstandard:
        warnings.append("non-standard residues present: " + ", ".join(nonstandard)
                        + " (may need custom parameters in pdb2gmx)")

    if free_sg and bonds:
        warnings.append(f"{len(free_sg)} cysteine(s) without an SG partner: "
                        f"pdb2gmx may need explicit SS bonds or protonation")

    print("\n[7] Warnings")
    if not warnings:
        print("    none")
    for w in warnings:
        print(f"    - {w}")

    print("\n[OK] Inspection finished (chain inventory saved as CSV in "
          f"{cfg.RESULT_DIR}). Update config.py, then run: "
          "python pipeline.py docking")
