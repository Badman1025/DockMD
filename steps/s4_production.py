import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import PROJECT_DIR, run, load_state, get_path
import config as cfg

GMX = cfg.TOOLS["gmx"]


def main():
    md_dir = get_path(load_state(), "md_dir", what="MD working directory")
    if not os.path.isdir(md_dir):
        sys.exit(f"!! Working directory not found: {md_dir}")

    dst_mdp = os.path.join(md_dir, "md.mdp")
    if not os.path.exists(dst_mdp):
        shutil.copy(os.path.join(PROJECT_DIR, "mdp", "md.mdp"), dst_mdp)
        print("[*] Copied md.mdp to working directory")
    else:
        print("[*] md.mdp already exists (skipped; delete it manually to update)")

    md_gro = os.path.join(md_dir, "md.gro")
    if os.path.exists(md_gro):
        print(f"\n[skip] {md_gro} already exists, skipping production MD.")
        print("       To rerun, delete the md.* files first.")
        return

    print("\n[1/2] grompp: generating production run input ...")
    run([GMX, "grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt",
         "-p", "topol.top", "-o", "md.tpr", "-maxwarn", "1"], cwd=md_dir)

    mdrun_cmd = [GMX, "mdrun", "-v", "-deffnm", "md"]
    md_cpt = os.path.join(md_dir, "md.cpt")
    if os.path.exists(md_cpt):
        print("[*] Found md.cpt checkpoint, resuming production MD ...")
        mdrun_cmd += ["-cpi", "md.cpt"]

    print("\n[2/2] mdrun: production MD (may take a long time)...")
    run(mdrun_cmd, cwd=md_dir)

    print("\n================ Production MD finished! ================")
    print(f"Trajectory: {os.path.join(md_dir, 'md.xtc')}")
    print("Next step: python pipeline.py analysis")
    print("Note: equilibration analyses use npt.xtc (see analysis --group dynamics/qc); "
          "for production trajectory analysis replace npt with md in the corresponding gmx commands.")


if __name__ == "__main__":
    main()
