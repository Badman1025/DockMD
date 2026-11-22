import glob
import os
import re
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config as cfg
from utils import PROJECT_DIR, run, save_state

PYMOL = cfg.TOOLS["pymol"]
OBABEL = cfg.TOOLS["obabel"]
VINA = cfg.TOOLS["vina"]


def parse_vina_mode1(log_file):
    try:
        txt = open(log_file).read()
        m = re.search(r'-----+\+------------.*?\n\s*1\s+(-?\d+\.\d+)', txt, re.S)
        if m:
            return float(m.group(1))
    except Exception:
        pass
    return None


def detect_pocket_center(pdb_path):
    het_xyz, heavy_xyz = [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                het_xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elif line.startswith("ATOM") and line[76:78].strip() != "H":
                heavy_xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if het_xyz:
        n = len(het_xyz)
        return (sum(p[0] for p in het_xyz)/n, sum(p[1] for p in het_xyz)/n, sum(p[2] for p in het_xyz)/n), "ligand"
    if heavy_xyz:
        n = len(heavy_xyz)
        return (sum(p[0] for p in heavy_xyz)/n, sum(p[1] for p in heavy_xyz)/n, sum(p[2] for p in heavy_xyz)/n), "protein"
    raise RuntimeError("No atoms found in PDB!")


def normalize_chains():
    chains = "".join(cfg.PROTEIN_CHAINS.split()).upper()
    if not chains:
        sys.exit("!! PROTEIN_CHAINS in config.py is empty (e.g. 'A' or 'AB')")
    return chains


def main():
    chains = normalize_chains()
    chain_sel = "+".join(chains)
    print(f"[0] Target chains: {'+'.join(chains)} ({len(chains)} chain(s))")

    raw_dir = os.path.join(PROJECT_DIR, cfg.RAW_DATA_DIR)
    protein_pdb = os.path.join(raw_dir, cfg.PROTEIN_PDB)
    if not os.path.isfile(protein_pdb):
        sys.exit(f"!! Protein file not found: {protein_pdb}")

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    temp_dir = os.path.join(PROJECT_DIR, cfg.RESULT_DIR, cfg.DOCKING_TEMP, run_id)
    result_dir = os.path.join(PROJECT_DIR, cfg.RESULT_DIR, cfg.DOCKING_DIR, run_id)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(result_dir, exist_ok=True)

    print(f"\n[1] Extracting chains {chain_sel} and getting pocket center...")
    temp_chains_het = os.path.join(temp_dir, f"{cfg.PROTEIN_NAME}_{chains}_het.pdb")
    run([PYMOL, "-cq", "-d",
         f"load {protein_pdb}, prot; remove solvent; "
         f"remove not chain {chain_sel}; save {temp_chains_het}, prot"])

    if cfg.FIXED_CENTER:
        CX, CY, CZ = cfg.FIXED_CENTER
        src = "fixed"
    else:
        (CX, CY, CZ), src = detect_pocket_center(temp_chains_het)
    print(f"  Box center: ({CX:.3f}, {CY:.3f}, {CZ:.3f}) [source: {src}]")

    receptor_pdb = os.path.join(result_dir, f"{cfg.PROTEIN_NAME}_receptor.pdb")
    with open(temp_chains_het) as f, open(receptor_pdb, "w") as g:
        for line in f:
            if not line.startswith("HETATM"):
                g.write(line)

    print("\n[2] Starting docking...")
    run([OBABEL, receptor_pdb, "-O", os.path.join(temp_dir, "rec.pdbqt"), "-xr"])
    sdf_files = sorted(glob.glob(os.path.join(raw_dir, "*.sdf")))
    if not sdf_files:
        sys.exit(f"!! No *.sdf ligand files found in {raw_dir}")

    scores = []
    for idx, sdf in enumerate(sdf_files, 1):
        lig = os.path.splitext(os.path.basename(sdf))[0]
        print(f"[{idx}/{len(sdf_files)}] {lig}")
        lig_pdbqt = os.path.join(temp_dir, f"{lig}.pdbqt")
        run([OBABEL, sdf, "-O", lig_pdbqt, "-h", "-p", "7.4"])

        out_docked = os.path.join(result_dir, f"{lig}_docked.pdbqt")
        log_file = os.path.join(result_dir, f"vina_{lig}.log")
        conf = (f"receptor = {os.path.join(temp_dir, 'rec.pdbqt')}\n"
                f"ligand = {lig_pdbqt}\n"
                f"center_x = {CX:.3f}\ncenter_y = {CY:.3f}\ncenter_z = {CZ:.3f}\n"
                f"size_x = {cfg.BOX_SIZE}\nsize_y = {cfg.BOX_SIZE}\nsize_z = {cfg.BOX_SIZE}\n"
                f"out = {out_docked}\n"
                f"exhaustiveness = {cfg.EXHAUSTIVENESS}\n"
                f"seed = {cfg.VINA_SEED}\ncpu = {cfg.VINA_CPU}\n")
        open(os.path.join(temp_dir, f"conf_{lig}.txt"), "w").write(conf)
        run([VINA, "--config", os.path.join(temp_dir, f"conf_{lig}.txt"), "--log", log_file])

        score = parse_vina_mode1(log_file)
        scores.append((lig, score))
        run([OBABEL, out_docked, "-O", os.path.join(result_dir, f"{lig}_docked.pdb"), "-f", "1", "-l", "1"])

    with open(os.path.join(result_dir, "scores.csv"), "w") as f:
        f.write("ligand,mode1_affinity\n")
        for name, s in scores:
            f.write(f"{name},{s}\n")

    save_state(docking_dir=result_dir, protein_name=cfg.PROTEIN_NAME,
               protein_chains=chains)
    print(f"\n[OK] Docking finished! Check {result_dir}/scores.csv and xxx_docked.pdb")


if __name__ == "__main__":
    main()
