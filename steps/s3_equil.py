import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import PROJECT_DIR, run, load_state, get_path
import config as cfg

GMX = cfg.TOOLS["gmx"]

MDP_DIR = os.path.join(PROJECT_DIR, "mdp")
MDP_FILES = ["ions.mdp", "em.mdp", "nvt.mdp", "npt.mdp"]


def main():
    md_dir = get_path(load_state(), "md_dir", what="MD working directory")
    if not os.path.isdir(md_dir):
        sys.exit(f"!! Working directory not found: {md_dir}")

    print(f"[*] Checking and copying MDP files to: {md_dir}")
    for fname in MDP_FILES:
        dst = os.path.join(md_dir, fname)
        if not os.path.exists(dst):
            shutil.copy(os.path.join(MDP_DIR, fname), dst)
            print(f"  + Copied {fname}")
        else:
            print(f"  - {fname} already exists (skipped; delete it manually to update)")

    print("\n================ Starting equilibration ================")

    ions_gro = os.path.join(md_dir, "complex_ions.gro")
    if os.path.exists(ions_gro):
        print(f"\n[1/4] {ions_gro} already exists, skipping ion addition (avoids duplicate ions)")
        print("     To redo: delete complex_ions.gro and ions.tpr, then restore topol.top")
    else:
        print("\n[1/4] Adding neutralizing ions (0.15M NaCl)...")
        run([GMX, "grompp", "-f", "ions.mdp", "-c", "complex_solv.gro",
             "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "1"], cwd=md_dir)
        run([GMX, "genion", "-s", "ions.tpr", "-o", "complex_ions.gro",
             "-p", "topol.top", "-pname", "NA", "-nname", "CL",
             "-neutral", "-conc", "0.15"], cwd=md_dir, stdin_str="SOL\n")

    print("\n[2/4] Energy minimization (EM)...")
    run([GMX, "grompp", "-f", "em.mdp", "-c", "complex_ions.gro",
         "-p", "topol.top", "-o", "em.tpr", "-maxwarn", "1"], cwd=md_dir)
    run([GMX, "mdrun", "-v", "-deffnm", "em"], cwd=md_dir)

    print("\n[3/4] NVT equilibration (100 ps, 300 K)...")
    run([GMX, "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro",
         "-p", "topol.top", "-o", "nvt.tpr", "-maxwarn", "1"], cwd=md_dir)
    run([GMX, "mdrun", "-v", "-deffnm", "nvt"], cwd=md_dir)

    print("\n[4/4] NPT equilibration (1 ns, 1 bar)...")
    run([GMX, "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro",
         "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr", "-maxwarn", "1"], cwd=md_dir)
    run([GMX, "mdrun", "-v", "-deffnm", "npt"], cwd=md_dir)

    print("\n================ Equilibration finished! ================")
    print(f"Please check npt.gro in {md_dir} for structural integrity.")
    print("Then run production MD: python pipeline.py md")


if __name__ == "__main__":
    main()
