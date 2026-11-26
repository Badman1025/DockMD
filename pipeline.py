import argparse
import json

import utils


def cmd_status(_args):
    state = utils.load_state()
    if not state:
        print("state.json does not exist. No steps have been run yet.")
        return
    print(json.dumps(state, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser(description="DockMD: docking -> MD pipeline")
    sub = parser.add_subparsers(dest="command", required=True)

    for name in ("docking", "prep", "equil", "md"):
        p = sub.add_parser(name)
        p.set_defaults(func=name)

    p_all = sub.add_parser("all", help="run docking -> prep -> equil in sequence")
    p_all.set_defaults(func="all")

    p_ana = sub.add_parser("analysis",
                           help="NPT equilibration-level analysis (run after equil)")
    p_ana.add_argument("--group", default="all",
                       choices=["all", "qc", "dynamics", "fel2d", "fel3d"])
    p_ana.set_defaults(func="analysis")

    p_mda = sub.add_parser("mdanalysis",
                           help="production MD-level analysis (run after md)")
    p_mda.add_argument("--group", default="all",
                       choices=["all", "dynamics", "fel2d", "fel3d"])
    p_mda.set_defaults(func="mdanalysis")

    sub.add_parser("status").set_defaults(func="status")

    p_info = sub.add_parser("info",
                            help="inspect chains, ligand and pocket assignment")
    p_info.add_argument("--pdb", default=None,
                        help="protein PDB (default: config PROTEIN_PDB)")
    p_info.add_argument("--ligand", default=None,
                        help="ligand file .sdf/.mol/.mol2/.pdb (default: auto-detect)")
    p_info.add_argument("--ligand-chain", default=None,
                        help="treat this chain of the PDB as the ligand (peptide/protein partner)")
    p_info.set_defaults(func="info")

    p_rep = sub.add_parser("report",
                           help="build report data + self-contained HTML report from analysis outputs")
    p_rep.add_argument("--level", default="all", choices=["all", "npt", "md"])
    p_rep.add_argument("--ai", action="store_true",
                       help="append LLM interpretation section (needs API key env)")
    p_rep.set_defaults(func="report")

    args = parser.parse_args()

    if args.func == "status":
        cmd_status(args)
        return

    if args.func == "info":
        from steps import s0_info
        s0_info.main(pdb_path=args.pdb, ligand_path=args.ligand,
                     ligand_chain=args.ligand_chain)
        return

    if args.func == "docking":
        from steps import s1_docking
        s1_docking.main()
    elif args.func == "prep":
        from steps import s2_prep
        s2_prep.main()
    elif args.func == "equil":
        from steps import s3_equil
        s3_equil.main()
    elif args.func == "md":
        from steps import s4_production
        s4_production.main()
    elif args.func == "analysis":
        from steps import s5_analysis
        s5_analysis.main(group=args.group, level="npt")
    elif args.func == "mdanalysis":
        from steps import s5_analysis
        s5_analysis.main(group=args.group, level="md")
    elif args.func == "report":
        from ai import report_data
        report_data.main(level=args.level, use_ai=args.ai)
    elif args.func == "all":
        from steps import s1_docking, s2_prep, s3_equil
        s1_docking.main()
        s2_prep.main()
        s3_equil.main()
        print("\n[all] docking -> prep -> equil finished.")
        print("      production MD takes a long time; run it manually when ready: python pipeline.py md")

    if args.func not in ("status", "report"):
        print("\n[OK] Current step finished. Check progress: python pipeline.py status")


if __name__ == "__main__":
    main()
