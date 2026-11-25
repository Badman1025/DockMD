<div align="center">
  <img src="DockMD_logo.png" alt="DockMD logo" width="140"/>
</div>

**DockMD**

A cross-platform (Windows / Linux) pipeline that chains molecular docking and  
GROMACS molecular dynamics into one driver: docking, system preparation,  
equilibration, production MD, and analysis/plotting.

## Project Layout

```
DockMD/
├── pipeline.py          Main entry point (all steps are run from here)
├── config.py            Single configuration source: parameters + tool paths
├── utils.py             Cross-platform command execution + state.json handoff
├── ai/
│   ├── report_data.py   Deterministic statistics + single-file HTML report
│   └── interpret.py     Optional LLM interpretation with numeric verification
├── steps/
│   ├── s1_docking.py    AutoDock Vina docking (PyMOL, Open Babel, Vina)
│   ├── s2_prep.py       MD system preparation (PDBFixer, pdb2gmx, acpype)
│   ├── s3_equil.py      Ion addition + EM + NVT + NPT equilibration
│   ├── s4_production.py Production MD run (md.mdp)
│   └── s5_analysis.py   Analysis and plotting (QC / RMSD / Rg / RMSF / H-bonds / FEL)
├── mdp/                 GROMACS mdp templates (ions, em, nvt, npt, md)
├── 1.Raw_data/          ALL input files go here (protein PDB, ligand SDFs)
└── 2.Result/            ALL outputs go here (docking, MD setup, plots, state.json)
```

Steps are chained through `2.Result/state.json`, which records the docking  
result directory, the MD working directory, and the protein/ligand names. Each  
step is resumable: if its output files already exist, that step is skipped.

## Installation

Create the conda environment (name: `DockMD`):

```
conda env create -f environment.yml
conda activate DockMD
```

`environment.yml` installs everything from conda (conda-forge + bioconda  
channels): Python 3.10, numpy, scipy, matplotlib, plotly, pymol-open-source,  
openbabel, openmm + pdbfixer (used by prep), gromacs, `autodock-vina`  
(AutoDock Vina 1.2.x, provides the `vina` command), plus the pip packages  
`kaleido` (fel3d plotting) and `acpype`.

Verify the toolchain after installation:

```
pymol -cq -d "print('pymol ok')"
obabel -V
vina --version
acpype --version
gmx --version
python -c "import numpy, scipy, matplotlib, plotly, openmm, pdbfixer; print('python deps ok')"
```

Platform notes for the environment:

- `gromacs` from conda-forge is available on Linux and macOS. On native  
  Windows it is not provided; either install GROMACS from the official  
  installer and put `gmx.exe` on PATH, or run the whole pipeline inside WSL2  
  (recommended for MD). A step-by-step Windows installation guide:  
  <https://structbiotools.com/how-to-install-gromacs/#windows>.  
  If tools are installed outside conda, set their  
  absolute paths in `config.py` -> `TOOLS`.
- If `autodock-vina` cannot be resolved on your platform, download the  
  standalone AutoDock Vina binary and point `TOOLS["vina"]` at it.
- On a cluster, it is also common to keep this conda env for docking/analysis  
  and load GROMACS via the cluster module system (`module load gromacs`),  
  adjusting `TOOLS["gmx"]` accordingly.

## Prerequisites

External tools must be available (on PATH, or set as absolute paths in  
`config.py` -> `TOOLS`). Installing the conda environment from  
`environment.yml` (see Installation) provides all of them except a native  
Windows GROMACS:

| Tool   | Used by    | Purpose                       |
| ------ | ---------- | ----------------------------- |
| pymol  | s1_docking | Extract the configured chains from the PDB    |
| obabel | s1, s2     | PDB/PDBQT/MOL2 conversion     |
| vina   | s1_docking | Molecular docking             |
| acpype | s2_prep    | GAFF2 ligand parameterization |
| gmx    | s2-s5      | GROMACS                       |

## Input Files

Place these in the `1.Raw_data/` directory before running step 1:

- `3HOR.pdb` (or change `PROTEIN_PDB` in `config.py`)
- One or more ligand files `*.sdf`; the ligand name is the file name without  
  the extension and must match `LIGAND_NAME` in `config.py` for the MD steps
- The chain IDs set in `PROTEIN_CHAINS` must exist in the PDB; PDB chain IDs  
  are single characters, so `"AB"` means chains A and B (not a chain named  
  "AB")

## Configuration

Edit `config.py`:

- `PROTEIN_PDB`, `PROTEIN_NAME`, `LIGAND_NAME`
- `PROTEIN_CHAINS`: chain IDs of the receptor to keep, written as one string  
  of single-character PDB chain IDs. `"A"` = monomer, `"AB"` = heterodimer  
  (chains A + B), `"ABC"` = trimer, and so on. Docking extracts only these  
  chains; prep builds the topology with one `Protein_chain_X` entry per chain  
  (generated automatically by pdb2gmx); the MD working directory is named  
  `<PROTEIN>_<CHAINS>_<LIGAND>` so different stoichiometries never mix
- Docking parameters: `BOX_SIZE`, `EXHAUSTIVENESS`, `VINA_SEED`, `VINA_CPU`,  
  `FIXED_CENTER` (set a `(x, y, z)` tuple to override pocket auto-detection)
- `TOOLS`: per-platform executable paths; if a tool is not on PATH, write its  
  absolute path (e.g. `r"C:\gromacs\bin\gmx.exe"`)

## Usage

Run all commands from the `DockMD/` directory.

```
python pipeline.py docking          Step 1: docking
python pipeline.py prep             Step 2: MD system preparation
python pipeline.py equil            Step 3: ions + EM + NVT + NPT
python pipeline.py md               Step 4: production MD
python pipeline.py analysis         Step 5a: NPT-level analysis (run after equil)
python pipeline.py mdanalysis       Step 5b: production-level analysis (run after md)
python pipeline.py analysis --group qc         Equilibration QC energy plots
python pipeline.py analysis --group dynamics   RMSD / Rg / RMSF / H-bonds
python pipeline.py analysis --group fel2d      PCA + 2D free energy landscape
python pipeline.py analysis --group fel3d      3D free energy landscape
python pipeline.py mdanalysis --group dynamics RMSD / Rg / RMSF / H-bonds (md.xtc)
python pipeline.py mdanalysis --group fel2d    PCA + 2D FEL (md.xtc)
python pipeline.py mdanalysis --group fel3d    3D FEL (md.xtc)
python pipeline.py all              Steps 1-3 in sequence (does NOT auto-run md)
python pipeline.py status           Show current progress from state.json
python pipeline.py report           Build 12/13 report files from existing xvg data
python pipeline.py report --level md   Report only the production level
python pipeline.py report --ai      Append the LLM interpretation section
```

Recommended order: run the steps one at a time and inspect the outputs before  
continuing. `all` stops after equilibration on purpose, since production MD  
takes hours to days.

## Outputs

All outputs live under `2.Result/`.

| Step       | Location (under `2.Result/`)      | Key outputs                                        |
| ---------- | --------------------------------- | -------------------------------------------------- |
| docking    | `Docking_Results/<timestamp>/`    | `scores.csv`, `*_docked.pdb`, `*_receptor.pdb`     |
| prep       | `MD_Setup/<PROT>_<CHAINS>_<LIG>/` | `topol.top`, `complex_solv.gro`, `final_check.pdb` |
| equil      | `MD_Setup/<PROT>_<CHAINS>_<LIG>/` | `em.*`, `nvt.*`, `npt.*`                           |
| md         | `MD_Setup/<PROT>_<CHAINS>_<LIG>/` | `md.xtc`, `md.gro`, `md.edr`                       |
| analysis   | `Analysis/<PROT>_<CHAINS>_<LIG>/` | files 01-06 in the table below (npt level)         |
| mdanalysis | `Analysis/<PROT>_<CHAINS>_<LIG>/` | files 07-11 in the table below (md level)          |

Analysis is split into two levels, run in this order:

1. `analysis` (NPT level) - available right after `equil`; uses `npt.tpr` /  
   `npt.xtc`. This is a stability quick-check on the 1 ns equilibration  
   trajectory, not a production result.
2. `mdanalysis` (MD level) - available after `md`; uses `md.tpr` / `md.xtc`.  
   This is the formal analysis of the full production trajectory.

All figures are collected in one folder per system (`cfg.ANALYSIS_DIR`,  
default `2.Result/Analysis/`) and distinguished by file name: files 01-06  
belong to the NPT level, files 07-11 (with `md` in the name) belong to the  
production level:

| File                            | Group    | Level | Content                                  |
| ------------------------------- | -------- | ----- | ---------------------------------------- |
| `01_Equilibration_QC.png`       | qc       | npt   | 2x2 EM/NVT/NPT energy, pressure, density |
| `02_Dynamics_Analysis.pdf`      | dynamics | npt   | RMSD, Rg, RMSF, protein-ligand H-bonds   |
| `03_FEL_Plot.png`               | fel2d    | npt   | 2D free energy landscape (standard)      |
| `04_FEL_Plot_Advanced.png`      | fel2d    | npt   | 2D FEL, Gaussian-smoothed                |
| `05_FEL_3D_Advanced.png`        | fel3d    | npt   | 3D free energy surface (static)          |
| `06_FEL_3D_Interactive.html`    | fel3d    | npt   | 3D interactive surface (Plotly)          |
| `07_md_Dynamics_Analysis.pdf`   | dynamics | md    | RMSD, Rg, RMSF, H-bonds on md.xtc        |
| `08_md_FEL_Plot.png`            | fel2d    | md    | 2D free energy landscape (standard)      |
| `09_md_FEL_Plot_Advanced.png`   | fel2d    | md    | 2D FEL, Gaussian-smoothed                |
| `10_md_FEL_3D_Advanced.png`     | fel3d    | md    | 3D free energy surface (static)          |
| `11_md_FEL_3D_Interactive.html` | fel3d    | md    | 3D interactive surface (Plotly)          |

Intermediate data (`*.xvg`, `<level>_center.xtc`, `index_<level>.ndx`) stays  
in the MD working directory; intermediate file names carry the level prefix  
so the two analysis levels never overwrite each other.

`python pipeline.py report` additionally writes two files into the same  
analysis folder (from existing xvg data only, no GROMACS re-run):

| File                     | Content                                             |
| ------------------------ | --------------------------------------------------- |
| `12_AI_Report_Data.json` | Deterministic statistics (RMSD/Rg/RMSF/H-bonds/FEL) |
| `13_AI_Report.html`      | English HTML report (inline SVG charts; header logo |
|                          | video embedded as base64; FEL media linked)         |

## Report Layer

`ai/report_data.py` parses the xvg files produced by `analysis`/`mdanalysis`  
and computes deterministic statistics: RMSD trend classification (plateau /  
slow drift / evident drift at 0.05 / 0.2 A per ns thresholds), Rg change,  
top flexible residues, H-bond occupancy, FEL basin count below 2 kJ/mol,  
and equilibration QC summaries.

```
python pipeline.py report
python pipeline.py report --ai   # optional LLM interpretation section
```

- `13_AI_Report.html` is a single-file English HTML report (light theme,  
  print-friendly, professional MD terminology) with stat cards, inline SVG  
  charts redrawn from the same xvg data, the outputs checklist, an  
  auto-generated interpretation section, and an optional AI section.
- FEL media are referenced as relative links per level when present: the  
  advanced 2D landscape PNG (04 / 09) and the interactive 3D FEL (06 / 11).  
  The report stays at ~6 MB because nothing beyond the logo video is  
  embedded; the links resolve as long as the report file is kept in the  
  analysis output directory.
- The interpretation section is written by the deterministic script at  
  report time (RMSD trend verdict, Rg stability, RMSF hotspots with a  
  rebuilt-loop caveat, H-bond occupancy classes, FEL basin summary, QC  
  pass/fail against the mdp target temperature). No AI image reading or  
  manual backfill is involved; the HTML is complete when generated.
- System-aware wording: QC targets (`ref_t`/`ref_p`) are parsed from  
  `mdout.mdp` instead of assuming 300 K; trajectories shorter than 10 ns  
  get a short-track disclaimer instead of a trend verdict; FEL results  
  below 5000 frames are marked as non-quantitative.
- Time axes are unit-aware: xvg headers declaring "Time (ps)" are converted  
  to ns in charts and statistics, so pre-fix files (gyrate without `-tu ns`)  
  and post-fix files are both handled correctly.
- `--ai` appends an LLM interpretation section (`ai/interpret.py`, stdlib  
  only). Provider, model and key are configured in `config.py`:

  ```python
  LLM_BASE_URL    = "https://api.siliconflow.cn/v1"
  LLM_MODEL       = "Qwen/Qwen2.5-7B-Instruct"
  LLM_API_KEY_ENV = "DASHSCOPE_API_KEY"
  ```

  `LLM_API_KEY_ENV` is the NAME of the environment variable holding the  
  API key; set it before generating the report:

  ```
  export DASHSCOPE_API_KEY=sk-...
  python pipeline.py report --ai
  ```

  The key is never written to any file in the repository. Any OpenAI-  
  compatible endpoint works — to switch providers, edit these three  
  lines only.  
  The model ID must match the provider's exact callable ID (e.g. the  
  marketplace ID `Qwen/Qwen2.5-7B-Instruct` on SiliconFlow). The model  
  receives the JSON statistics (not images), must cite only numbers  
  present in the JSON, and every number it writes is programmatically  
  verified against the JSON — sentences with unverifiable numbers are  
  dropped. Failure modes (missing key, network error, all-sentences-  
  rejected) degrade silently to the deterministic report; the run never  
  fails because of the AI layer.
- All interpretation text (deterministic and AI-assisted) is descriptive and  
  requires human review before use in formal materials. The report layer  
  never touches the physics core: no mdp editing, no parameter changes.

## Running in the Background

Production MD (step 4) takes hours to days. Run it detached with `nohup`  
from a Linux/WSL2 terminal:

```
cd DockMD
conda activate DockMD
nohup python -u pipeline.py md > md_pipeline.log 2>&1 &
tail -f md_pipeline.log
```

- `-u` disables Python output buffering so the log file updates in real time;
- the pipeline can be started from any directory, but keep `md_pipeline.log`  
  where you like; `Ctrl+C` / closing the terminal will NOT kill the job once  
  started with `nohup ... &`;
- check GROMACS progress in the run log:  
  `tail -f 2.Result/MD_Setup/<PROT>_<CHAINS>_<LIG>/md.log`  
  (look for the step counter and "Performance" lines);
- stop a running job: `pkill -f "gmx mdrun"` then `pkill -f "pipeline.py md"`.

If the job is interrupted, simply re-run `python pipeline.py md`: it detects  
`md.cpt` and resumes from the last checkpoint with `mdrun -cpi md.cpt`  
(no restart from scratch). Only after the run completes (`md.gro` exists)  
does the step skip itself.

On native Windows there is no `nohup`; use WSL2 for long MD runs, or run the  
command in a terminal session you keep open.

## Resume Behavior

- `equil`: skips ion addition if `complex_ions.gro` exists (delete it, plus  
  `ions.tpr`, and restore `topol.top` to redo).

- `md`: skips if `md.gro` exists (delete the `md.*` files to rerun).
- `prep`: reads the docking directory from `state.json`; if absent it falls  
  back to the newest directory under `Docking_Results/`.

## Platform Notes

- The code is identical on Windows and Linux; all external commands are  
  launched as argument lists without a shell.
- For long MD runs (NVT/NPT/production), Linux or WSL2 is recommended:  
  native GROMACS on Windows has weaker GPU and MPI support. The same code and  
  `state.json` can be copied between machines; just adjust `TOOLS` paths.
- `analysis` (NPT level) only needs `npt.tpr` / `npt.xtc` from equilibration,  
  and `mdanalysis` (MD level) only needs `md.tpr` / `md.xtc` from the  
  production run, so plotting can be done on a different machine than the MD  
  run.
