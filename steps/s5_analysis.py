import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils import PROJECT_DIR, run, load_state, get_path
import config as cfg

GMX = cfg.TOOLS["gmx"]


def _setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    return plt


def group_qc(md_dir, level="npt"):
    def get_energy(xvg_file, term, out_file):
        print(f"Extracting data: {term} -> {out_file}")
        run([GMX, "energy", "-f", xvg_file, "-o", out_file],
            cwd=md_dir, stdin_str=f"{term}\n0\n", quiet=True)

    def plot_xvg(plt, xvg_file, title, ylabel, color="b", ax=None):
        data = []
        with open(os.path.join(md_dir, xvg_file), 'r') as f:
            for line in f:
                if line.startswith(('#', '@')):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    data.append((float(parts[0]), float(parts[1])))

        if not data:
            print(f"Warning: no data in {xvg_file}!")
            return

        x, y = zip(*data)

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5))

        ax.plot(x, y, color=color, linewidth=1.5)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Time (ps) / Step", fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.6)

    plt = _setup_matplotlib()

    print("Extracting data and plotting...")
    get_energy("em.edr", "Potential", "em_potential.xvg")
    get_energy("nvt.edr", "Temperature", "nvt_temp.xvg")
    get_energy("npt.edr", "Pressure", "npt_pres.xvg")
    get_energy("npt.edr", "Density", "npt_dens.xvg")
    get_energy("npt.edr", "Temperature", "npt_temp.xvg")

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('GROMACS Equilibration QC Plots', fontsize=16)

    plot_xvg(plt, "em_potential.xvg", "Energy Minimization",
             "Potential Energy (kJ/mol)", color="purple", ax=axs[0, 0])
    axs[0, 0].set_xlabel("Step")
    plot_xvg(plt, "nvt_temp.xvg", "NVT Equilibration",
             "Temperature (K)", color="red", ax=axs[0, 1])
    plot_xvg(plt, "npt_pres.xvg", "NPT Equilibration",
             "Pressure (bar)", color="green", ax=axs[1, 0])
    plot_xvg(plt, "npt_dens.xvg", "NPT Equilibration",
             "Density (kg/m^3)", color="blue", ax=axs[1, 1])

    plt.tight_layout()
    out_png = os.path.join(_analysis_dir(md_dir), "01_Equilibration_QC.png")
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"\nPlotting finished! Saved to: {out_png}")


def _analysis_dir(md_dir):
    out = os.path.join(PROJECT_DIR, cfg.RESULT_DIR, cfg.ANALYSIS_DIR,
                       os.path.basename(md_dir))
    os.makedirs(out, exist_ok=True)
    return out


def _resname_from_itp(md_dir):
    for name in (f"{cfg.LIGAND_NAME}_clean.itp", f"{cfg.LIGAND_NAME}_GMX.itp"):
        path = os.path.join(md_dir, f"{cfg.LIGAND_NAME}.acpype", name)
        if os.path.isfile(path):
            with open(path) as f:
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
    return None


def _lig_resname(md_dir):
    r = _resname_from_itp(md_dir)
    if r:
        return r
    state = load_state()
    if state.get("lig_resname"):
        return state["lig_resname"]
    topol = os.path.join(md_dir, "topol.top")
    if os.path.isfile(topol):
        with open(topol) as f:
            in_mols = False
            for line in f:
                s = line.strip()
                if s.startswith("[ molecules ]"):
                    in_mols = True
                    continue
                if in_mols:
                    if s.startswith("["):
                        break
                    if s and not s.startswith(";") and not s.startswith("#"):
                        name = s.split()[0]
                        if not name.startswith("Protein"):
                            return name
    return "UNK"


def _read_groups(ndx_path):
    groups = []
    with open(ndx_path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("[") and s.endswith("]"):
                groups.append(s[1:-1].strip())
    return groups


def _find_group_number(ndx_path, candidates):
    groups = _read_groups(ndx_path)
    for c in candidates:
        if c in groups:
            return groups.index(c)
    return None


def _pick_hbond_tool():
    probe = subprocess.run([GMX, "hbond-legacy", "-h"],
                           capture_output=True, text=True)
    if probe.returncode == 0:
        return "hbond-legacy", True
    return "hbond", False


def group_dynamics(md_dir, level="npt"):
    plt = _setup_matplotlib()

    resname = _lig_resname(md_dir)
    tpr = f"{level}.tpr"
    traj = f"{level}.xtc"
    center = f"{level}_center.xtc"
    print(f"\n[1/2] Running GROMACS trajectory analysis on {traj} "
          f"(ligand residue: {resname}) ...")
    run([GMX, "rms", "-s", tpr, "-f", traj, "-o", f"{level}_rmsd.xvg", "-tu", "ns"],
        cwd=md_dir, stdin_str="Backbone\nProtein\n", quiet=True)
    run([GMX, "gyrate", "-s", tpr, "-f", traj, "-o", f"{level}_gyrate.xvg", "-tu", "ns"],
        cwd=md_dir, stdin_str="Protein\n", quiet=True)
    run([GMX, "trjconv", "-s", tpr, "-f", traj,
         "-o", center, "-pbc", "mol", "-center"],
        cwd=md_dir, stdin_str="Protein\nSystem\n", quiet=True)
    run([GMX, "rmsf", "-s", tpr, "-f", center, "-o", f"{level}_rmsf.xvg", "-res"],
        cwd=md_dir, stdin_str="Protein\n", quiet=True)
    run([GMX, "make_ndx", "-f", tpr, "-o", f"index_{level}.ndx"],
        cwd=md_dir, stdin_str="q\n", quiet=True)
    lig_num = _find_group_number(os.path.join(md_dir, f"index_{level}.ndx"),
                                 [f"r_{resname}", resname])
    if lig_num is None:
        groups = _read_groups(os.path.join(md_dir, f"index_{level}.ndx"))
        sys.exit(f"!! Ligand group for residue '{resname}' not found in "
                 f"index_{level}.ndx. Available groups: {groups}")
    hb_tool, is_legacy = _pick_hbond_tool()
    if is_legacy:
        print(f"Using gmx {hb_tool} (legacy hbond, same as original script)")
        run([GMX, hb_tool, "-s", tpr, "-f", traj, "-n", f"index_{level}.ndx",
             "-num", f"{level}_hbnum.xvg", "-tu", "ns", "-r", "0.35"],
            cwd=md_dir, stdin_str=f"1\n{lig_num}\n", quiet=True)
    else:
        print("Warning: gmx hbond-legacy not found in this GROMACS build; "
              "falling back to the rewritten gmx hbond. If the ligand has no "
              "donors/acceptors detected, install a GROMACS build that ships "
              "hbond-legacy (conda-forge builds include it).")
        run([GMX, hb_tool, "-s", tpr, "-f", traj, "-n", f"index_{level}.ndx",
             "-num", f"{level}_hbnum.xvg", "-tu", "ns"],
            cwd=md_dir, stdin_str=f"1\n{lig_num}\n", quiet=True)

    if level == "npt":
        pdf_path = os.path.join(_analysis_dir(md_dir), "02_Dynamics_Analysis.pdf")
    else:
        pdf_path = os.path.join(_analysis_dir(md_dir), "07_md_Dynamics_Analysis.pdf")
    print(f"\n[2/2] Plotting {os.path.basename(pdf_path)} ...")

    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages

    def read_xvg(filename):
        if not os.path.exists(os.path.join(md_dir, filename)):
            return None, None
        try:
            data = np.loadtxt(os.path.join(md_dir, filename), comments=['#', '@'])
            return data[:, 0], data[:, 1]
        except Exception:
            return None, None

    pdf = PdfPages(pdf_path)

    time_rmsd, rmsd = read_xvg(f"{level}_rmsd.xvg")
    if time_rmsd is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(time_rmsd, rmsd * 10, color='#1f77b4', linewidth=1.5)
        avg_rmsd = np.mean(rmsd[len(rmsd)//2:] * 10)
        ax.axhline(y=avg_rmsd, color='red', linestyle='--', linewidth=1.2, label=f'Avg: {avg_rmsd:.2f} A')
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Backbone RMSD (A)', fontsize=12)
        ax.set_title('Protein Backbone RMSD', fontsize=14, fontweight='bold')
        ax.legend(frameon=False, fontsize=10)
        ax.grid(True, linestyle='--', alpha=0.6)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close()

    time_rg, rg = read_xvg(f"{level}_gyrate.xvg")
    if time_rg is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(time_rg, rg, color='#2ca02c', linewidth=1.5)
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Radius of Gyration (nm)', fontsize=12)
        ax.set_title('Protein Radius of Gyration (Rg)', fontsize=14, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.6)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close()

    resi_rmsf, rmsf = read_xvg(f"{level}_rmsf.xvg")
    if resi_rmsf is not None:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(resi_rmsf, rmsf * 10, color='#d62728', linewidth=1.2)
        ax.set_xlabel('Residue Number', fontsize=12)
        ax.set_ylabel('RMSF (A)', fontsize=12)
        ax.set_title('Residue Fluctuation (RMSF)', fontsize=14, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.6)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close()

    time_hb, hbonds = read_xvg(f"{level}_hbnum.xvg")
    if time_hb is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(time_hb, hbonds, color='#9467bd', linewidth=1.5)
        if np.sum(hbonds) == 0:
            ax.text(0.5, 0.5, 'No H-bonds found', transform=ax.transAxes, ha='center', va='center', color='red')
        ax.set_xlabel('Time (ns)', fontsize=12)
        ax.set_ylabel('Number of H-bonds', fontsize=12)
        ax.set_title('Protein-Ligand Hydrogen Bonds', fontsize=14, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.6)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close()

    pdf.close()
    print(f"Generated: {pdf_path}")


def group_fel2d(md_dir, level="npt"):
    plt = _setup_matplotlib()

    tpr = f"{level}.tpr"
    center = f"{level}_center.xtc"

    print("\n[1/2] Running GROMACS PCA analysis ...")
    run([GMX, "covar", "-s", tpr, "-f", center,
         "-o", f"{level}_eigenvalues.xvg", "-v", f"{level}_eigenvectors.trr",
         "-av", f"{level}_average.pdb", "-l", f"{level}_covar.log"],
        cwd=md_dir, stdin_str="Backbone\nProtein\n", quiet=True)
    run([GMX, "anaeig", "-s", tpr, "-f", center,
         "-v", f"{level}_eigenvectors.trr", "-2d", f"{level}_proj_2d.xvg",
         "-first", "1", "-last", "2"],
        cwd=md_dir, stdin_str="Backbone\nProtein\n", quiet=True)

    print("\n[2/2] Plotting 2D free energy landscape ...")

    import numpy as np
    from scipy.ndimage import gaussian_filter

    def read_2d_xvg(filename):
        path = os.path.join(md_dir, filename)
        if not os.path.exists(path):
            return None, None, None
        try:
            data = np.loadtxt(path, comments=['#', '@'])
            if data.shape[1] == 3:
                return data[:, 0], data[:, 1], data[:, 2]
            else:
                return None, data[:, 0], data[:, 1]
        except Exception:
            return None, None, None

    time, pc1, pc2 = read_2d_xvg(f"{level}_proj_2d.xvg")

    if pc1 is None:
        print(f"{level}_proj_2d.xvg not found or format not supported. "
              f"(Hint: run the dynamics group at the same level first to "
              f"generate {center})")
        return

    hist1, xedges1, yedges1 = np.histogram2d(pc1, pc2, bins=50, density=True)
    fe1 = np.zeros_like(hist1)
    valid1 = hist1 > 0
    fe1[valid1] = -np.log(hist1[valid1])
    fe1 = fe1 - np.min(fe1)
    fe1[~valid1] = np.nan

    fig1, ax1 = plt.subplots(figsize=(6, 5))
    im1 = ax1.imshow(fe1.T, origin='lower', aspect='auto',
                     extent=[xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]],
                     cmap='jet', interpolation='bilinear', vmax=8)
    cbar1 = fig1.colorbar(im1, ax=ax1, pad=0.03)
    cbar1.set_label('G (kJ/mol)', fontsize=12)
    ax1.set_xlabel('PC1', fontsize=12)
    ax1.set_ylabel('PC2', fontsize=12)
    ax1.set_title('Gibbs Energy Landscape', fontsize=14, fontweight='bold')
    fig1.tight_layout()
    if level == "npt":
        output_path1 = os.path.join(_analysis_dir(md_dir), "03_FEL_Plot.png")
        output_path2 = os.path.join(_analysis_dir(md_dir), "04_FEL_Plot_Advanced.png")
    else:
        output_path1 = os.path.join(_analysis_dir(md_dir), "08_md_FEL_Plot.png")
        output_path2 = os.path.join(_analysis_dir(md_dir), "09_md_FEL_Plot_Advanced.png")
    plt.savefig(output_path1, dpi=300)
    plt.close(fig1)
    print(f"Generated standard version: {output_path1}")

    hist2, xedges2, yedges2 = np.histogram2d(pc1, pc2, bins=50, density=True)
    fe2 = np.zeros_like(hist2)
    valid2 = hist2 > 0
    fe2[valid2] = -np.log(hist2[valid2])
    fe2 = fe2 - np.min(fe2)

    vmax_val = 8.0
    fe2[~valid2] = vmax_val
    fe2 = gaussian_filter(fe2, sigma=1.0)
    fe2 = np.clip(fe2, 0, vmax_val)

    fig2, ax2 = plt.subplots(figsize=(6.5, 5.5))
    im2 = ax2.imshow(fe2.T, origin='lower', aspect='auto',
                     extent=[xedges2[0], xedges2[-1], yedges2[0], yedges2[-1]],
                     cmap='jet', interpolation='bicubic', vmin=0, vmax=vmax_val)
    cbar2 = fig2.colorbar(im2, ax=ax2, pad=0.03)
    cbar2.set_label('Gibbs Free Energy (kJ/mol)', fontsize=13, labelpad=10)
    cbar2.ax.tick_params(labelsize=11)

    ax2.set_xlabel('PC1', fontsize=14, fontweight='bold')
    ax2.set_ylabel('PC2', fontsize=14, fontweight='bold')
    ax2.set_title('Gibbs Energy Landscape', fontsize=16, fontweight='bold', pad=15)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    fig2.tight_layout()
    plt.savefig(output_path2, dpi=300)
    plt.close(fig2)
    print(f"Generated advanced version: {output_path2}")


def group_fel3d(md_dir, level="npt"):
    import numpy as np
    import plotly.graph_objects as go
    from scipy.ndimage import gaussian_filter

    filepath = os.path.join(md_dir, f"{level}_proj_2d.xvg")

    if not os.path.exists(filepath):
        print(f"File not found: {filepath}. Run the fel2d group at the same "
              f"level first to generate the PCA projection.")
        return

    data = np.loadtxt(filepath, comments=['#', '@'])
    pc1 = data[:, 0]
    pc2 = data[:, 1]

    hist, xedges, yedges = np.histogram2d(pc1, pc2, bins=50, density=True)
    free_energy = np.zeros_like(hist)
    valid_bins = hist > 0
    free_energy[valid_bins] = -np.log(hist[valid_bins])

    vmax_val = 8.0
    free_energy[~valid_bins] = vmax_val
    free_energy = gaussian_filter(free_energy, sigma=1.0)
    free_energy = np.clip(free_energy, 0, vmax_val)

    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(x_centers, y_centers)
    Z = free_energy.T

    fig = go.Figure(data=[
        go.Surface(
            z=Z, x=X, y=Y,
            colorscale='Jet',
            cmin=0, cmax=vmax_val,
            colorbar=dict(title=dict(text='G (kJ/mol)', side='right'), len=0.7),
            name='Energy Surface',
            showscale=True
        ),
        go.Surface(
            z=np.full_like(Z, -0.5),
            x=X, y=Y,
            surfacecolor=Z,
            colorscale='Jet',
            cmin=0, cmax=vmax_val,
            showscale=False,
            opacity=0.8,
            contours=dict(
                z=dict(show=True, usecolormap=True, highlightcolor="black", project=dict(z=True))
            )
        )
    ])

    fig.update_layout(
        title=dict(
            text='3D Gibbs Free Energy Landscape',
            font=dict(size=22, family='Arial', weight='bold'),
            x=0.5
        ),
        scene=dict(
            xaxis_title='PC1',
            yaxis_title='PC2',
            zaxis_title='G (kJ/mol)',
            xaxis=dict(backgroundcolor="white", gridcolor="lightgrey", showbackground=True),
            yaxis=dict(backgroundcolor="white", gridcolor="lightgrey", showbackground=True),
            zaxis=dict(backgroundcolor="white", gridcolor="lightgrey", showbackground=True),
            camera=dict(eye=dict(x=1.5, y=-1.5, z=1.0))
        ),
        width=1000,
        height=850,
        margin=dict(l=50, r=50, b=50, t=80),
        template='plotly_white'
    )

    if level == "npt":
        output_png = os.path.join(_analysis_dir(md_dir), "05_FEL_3D_Advanced.png")
        output_html = os.path.join(_analysis_dir(md_dir), "06_FEL_3D_Interactive.html")
    else:
        output_png = os.path.join(_analysis_dir(md_dir), "10_md_FEL_3D_Advanced.png")
        output_html = os.path.join(_analysis_dir(md_dir), "11_md_FEL_3D_Interactive.html")

    fig.write_image(output_png, scale=2)
    fig.write_html(output_html)

    print(f"3D high-resolution static image generated: {output_png}")
    print(f"3D interactive HTML generated: {output_html}")


NPT_GROUPS = {
    "qc": ("Equilibration QC energy plots", group_qc),
    "dynamics": ("RMSD/Rg/RMSF/H-bonds", group_dynamics),
    "fel2d": ("2D free energy landscape", group_fel2d),
    "fel3d": ("3D free energy landscape", group_fel3d),
}

MD_GROUPS = {
    "dynamics": ("RMSD/Rg/RMSF/H-bonds", group_dynamics),
    "fel2d": ("2D free energy landscape", group_fel2d),
    "fel3d": ("3D free energy landscape", group_fel3d),
}


def main(group="all", level="npt"):
    md_dir = get_path(load_state(), "md_dir", what="MD working directory")
    if not os.path.isdir(md_dir):
        sys.exit(f"!! Working directory not found: {md_dir}")

    groups = NPT_GROUPS if level == "npt" else MD_GROUPS
    if group == "all":
        todo = list(groups)
    elif group in groups:
        todo = [group]
    else:
        sys.exit(f"!! Group '{group}' is not available at level '{level}'. "
                 f"Available groups: {list(groups)}")
    if level == "md" and not os.path.isfile(os.path.join(md_dir, "md.xtc")):
        sys.exit("!! md.xtc not found in the MD working directory. "
                 "Run 'python pipeline.py md' first.")
    if level == "npt" and not os.path.isfile(os.path.join(md_dir, "npt.xtc")):
        sys.exit("!! npt.xtc not found in the MD working directory. "
                 "Run 'python pipeline.py equil' first.")

    for g in todo:
        title, fn = groups[g]
        print(f"\n================ Analysis [{level}]: {title} ({g}) ================")
        fn(md_dir, level)

    print(f"\n[OK] Analysis finished. All figures are in "
          f"{os.path.join(PROJECT_DIR, cfg.RESULT_DIR, cfg.ANALYSIS_DIR, os.path.basename(md_dir))}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="all",
                    choices=["all", "qc", "dynamics", "fel2d", "fel3d"])
    ap.add_argument("--level", default="npt", choices=["npt", "md"])
    ns = ap.parse_args()
    main(group=ns.group, level=ns.level)
