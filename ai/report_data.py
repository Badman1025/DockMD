import base64
import html as _html
import json
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

from utils import PROJECT_DIR, load_state
import config as cfg

RMSD_PLATEAU = 0.05
RMSD_DRIFT = 0.2
FEL_VMAX = 8.0
FEL_BASIN_CUTOFF = 2.0
FEL_MIN_FRAMES = 5000
TOP_FLEXIBLE = 8
LOGO_VIDEO = os.path.join(PROJECT_DIR, "Local_PyMOL_views", "dockVD_logo.mp4")

ANALYSIS_ROOT = os.path.join(PROJECT_DIR, cfg.RESULT_DIR, cfg.ANALYSIS_DIR)

EXPECTED_OUTPUTS = {
    "npt": [
        "01_Equilibration_QC.png",
        "02_Dynamics_Analysis.pdf",
        "03_FEL_Plot.png",
        "04_FEL_Plot_Advanced.png",
        "05_FEL_3D_Advanced.png",
        "06_FEL_3D_Interactive.html",
    ],
    "md": [
        "07_md_Dynamics_Analysis.pdf",
        "08_md_FEL_Plot.png",
        "09_md_FEL_Plot_Advanced.png",
        "10_md_FEL_3D_Advanced.png",
        "11_md_FEL_3D_Interactive.html",
    ],
}

CSS = """
:root { --ink:#1f2937; --muted:#6b7280; --line:#e5e7eb; --accent:#1d4ed8; }
* { box-sizing: border-box; }
body { font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
       color: var(--ink); background: #f1f5f9; margin: 0; }
.wrap { max-width: 920px; margin: 24px auto 64px; background: #fff;
        border: 1px solid var(--line); border-radius: 12px; padding: 32px 40px 48px; }
header { display: flex; align-items: center; gap: 24px; border-bottom: 3px solid var(--accent);
         padding-bottom: 20px; margin-bottom: 8px; }
header .title { flex: 1; }
h1 { font-size: 24px; margin: 0 0 6px; letter-spacing: 0.5px; }
.meta { font-size: 12.5px; color: var(--muted); line-height: 1.8; }
.meta b { color: var(--ink); font-weight: 600; }
video.logo { width: 220px; border-radius: 10px; background: #000; flex-shrink: 0; }
.notice { background: #eff6ff; border: 1px solid #bfdbfe; color: #1e40af; border-radius: 8px;
          padding: 10px 14px; font-size: 13px; margin: 16px 0 8px; line-height: 1.7; }
h2 { font-size: 17px; margin: 34px 0 12px; padding-left: 10px;
     border-left: 4px solid var(--accent); line-height: 1.4; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; }
.card { background: #f8fafc; border: 1px solid var(--line); border-radius: 8px; padding: 10px 14px; }
.card-label { font-size: 11.5px; color: var(--muted); margin-bottom: 4px; }
.card-value { font-size: 16px; font-weight: 600; }
.card-value .sub { font-size: 12px; font-weight: 400; color: var(--muted); margin-left: 6px; }
.charts { display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 14px; margin-top: 14px; }
figure { margin: 0; border: 1px solid var(--line); border-radius: 8px; padding: 10px 10px 4px; }
figure img { width: 100%; height: auto; display: block; }
.file-links { list-style: none; margin: 10px 0 0; padding: 0; }
.file-links li { margin: 6px 0; font-size: 13.5px; }
.file-links a { color: var(--accent); font-weight: 600; text-decoration: none; }
.file-links a:hover { text-decoration: underline; }
.file-links .hint { color: var(--muted); font-size: 12px; }
figcaption { font-size: 12px; color: var(--muted); text-align: center; padding: 4px 0 6px; }
ul.plain { margin: 6px 0; padding-left: 20px; line-height: 1.9; font-size: 14px; }
.interp { border: 1px solid var(--line); border-radius: 8px; padding: 14px 18px; margin-bottom: 12px; }
.interp h3 { font-size: 14px; margin: 0 0 6px; }
.interp .src { font-size: 11.5px; color: var(--muted); margin-bottom: 6px; }
.interp-body { font-size: 14px; line-height: 1.9; }
.interp-body b { color: #15803d; }
.notice-amber { background: #fffbeb; border: 1px solid #fde68a; color: #92400e; border-radius: 8px;
                padding: 10px 14px; font-size: 13px; margin: 16px 0 8px; line-height: 1.7; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { border: 1px solid var(--line); padding: 5px 10px; text-align: left; }
th { background: #f8fafc; }
.ok { color: #15803d; } .miss { color: #b45309; }
@media print { body { background: #fff; } .wrap { border: none; margin: 0; } }
"""


def _read_xvg(path):
    if not os.path.isfile(path):
        return None
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(("#", "@")):
                continue
            try:
                rows.append([float(p) for p in s.split()])
            except ValueError:
                continue
    return rows or None


def _series(rows, scale=1.0):
    x = np.array([r[0] for r in rows], dtype=float)
    y = np.array([r[1] for r in rows], dtype=float) * scale
    return x, y


def _seg_stats(x, y, lo, hi):
    n = len(y)
    a = int(n * lo)
    b = max(int(n * hi), a + 1)
    return float(np.mean(y[a:b])), float(np.std(y[a:b]))


def _time_axis(md_dir, level, fname):
    rows = _read_xvg(os.path.join(md_dir, fname))
    if not rows:
        return None, None
    t = np.array([r[0] for r in rows], dtype=float)
    if _xvg_xunit(os.path.join(md_dir, fname)) == "ps":
        t = t / 1000.0
    return t, rows


def _mdp_targets(md_dir):
    targets = {"temperature": 300.0, "pressure": 1.0, "source": "defaults (mdout.mdp not parsed)"}
    path = os.path.join(md_dir, "mdout.mdp")
    if not os.path.isfile(path):
        return targets
    try:
        with open(path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.split(";")[0].split("#")[0].strip()
                if "=" not in s:
                    continue
                k, v = s.split("=", 1)
                k = k.strip().lower()
                parts = v.split()
                if not parts:
                    continue
                if k == "ref_t":
                    targets["temperature"] = float(parts[0])
                elif k == "ref_p":
                    targets["pressure"] = float(parts[0])
        targets["source"] = "mdout.mdp"
    except (OSError, ValueError):
        pass
    return targets


def _rmsd_block(md_dir, level):
    t, rows = _time_axis(md_dir, level, f"{level}_rmsd.xvg")
    if t is None:
        return None
    y = np.array([r[1] for r in rows], dtype=float) * 10.0
    n = len(y)
    init, _ = _seg_stats(t, y, 0.0, 0.05)
    late, late_sd = _seg_stats(t, y, 0.8, 1.0)
    half = n // 2
    slope = float(np.polyfit(t[half:], y[half:], 1)[0])
    imax = int(np.argmax(y))
    if abs(slope) < RMSD_PLATEAU:
        trend = "plateau"
    elif abs(slope) < RMSD_DRIFT:
        trend = "slow_drift"
    else:
        trend = "evident_drift"
    span = float(t[-1] - t[0])
    return {
        "unit": "Angstrom",
        "time_unit": "ns",
        "frames": n,
        "time_span_ns": round(span, 3),
        "short_track": bool(span < 10.0),
        "initial_mean": round(init, 3),
        "late_mean": round(late, 3),
        "late_std": round(late_sd, 3),
        "late_window_ns": [round(float(t[int(n * 0.8)]), 2), round(float(t[-1]), 2)],
        "max": round(float(y[imax]), 3),
        "max_time_ns": round(float(t[imax]), 2),
        "tail_slope_A_per_ns": round(slope, 4),
        "trend_class": trend,
    }


def _gyrate_block(md_dir, level):
    rows = _read_xvg(os.path.join(md_dir, f"{level}_gyrate.xvg"))
    if not rows:
        return None
    t, y = _series(rows)
    early, _ = _seg_stats(t, y, 0.0, 0.1)
    late, late_sd = _seg_stats(t, y, 0.8, 1.0)
    change = (late - early) / early * 100.0 if early else 0.0
    return {
        "unit": "nm",
        "frames": len(y),
        "early_mean": round(early, 3),
        "late_mean": round(late, 3),
        "late_std": round(late_sd, 3),
        "change_percent": round(change, 2),
    }


def _rmsf_block(md_dir, level):
    rows = _read_xvg(os.path.join(md_dir, f"{level}_rmsf.xvg"))
    if not rows:
        return None
    res = np.array([r[0] for r in rows], dtype=float)
    val = np.array([r[1] for r in rows], dtype=float) * 10.0
    order = np.argsort(val)[::-1]
    top = [{"residue": int(res[i]), "rmsf_A": round(float(val[i]), 2)} for i in order[:TOP_FLEXIBLE]]
    mean = float(np.mean(val))
    sd = float(np.std(val))
    return {
        "unit": "Angstrom",
        "n_residues": len(val),
        "mean": round(mean, 3),
        "std": round(sd, 3),
        "max": round(float(val.max()), 3),
        "top_flexible_residues": top,
        "n_above_mean_plus_2sd": int((val > mean + 2 * sd).sum()),
        "note": "residue numbers are GROMACS selection indices and may be offset from PDB numbering",
    }


def _hbnum_block(md_dir, level):
    t, rows = _time_axis(md_dir, level, f"{level}_hbnum.xvg")
    if t is None:
        return None
    y = np.array([r[1] for r in rows], dtype=float)
    y_int = np.round(y).astype(int)
    return {
        "unit": "count",
        "time_unit": "ns",
        "frames": len(y),
        "time_span_ns": round(float(t[-1] - t[0]), 3),
        "mean": round(float(np.mean(y)), 2),
        "max": int(y_int.max()),
        "fraction_frames_with_at_least_one": round(float(np.mean(y_int >= 1)), 4),
    }


def _fel_block(md_dir, level):
    rows = _read_xvg(os.path.join(md_dir, f"{level}_proj_2d.xvg"))
    if not rows or len(rows[0]) < 2:
        return None
    try:
        from scipy.ndimage import gaussian_filter, label
    except ImportError:
        gaussian_filter = None
        label = None
    pc1 = np.array([r[0] for r in rows], dtype=float)
    pc2 = np.array([r[1] for r in rows], dtype=float)
    hist, xe, ye = np.histogram2d(pc1, pc2, bins=50)
    bin_area = ((xe[-1] - xe[0]) / 50.0) * ((ye[-1] - ye[0]) / 50.0)
    total = hist.sum()
    if total > 0 and bin_area > 0:
        hist = hist / (total * bin_area)
    fe = np.zeros_like(hist)
    valid = hist > 0
    fe[valid] = -np.log(hist[valid])
    fe[valid] -= fe[valid].min()
    fe[~valid] = FEL_VMAX
    if gaussian_filter is not None:
        fe = gaussian_filter(fe, sigma=1.0)
    fe = np.clip(fe, 0, FEL_VMAX)
    xc = (xe[:-1] + xe[1:]) / 2.0
    yc = (ye[:-1] + ye[1:]) / 2.0
    X, Y = np.meshgrid(xc, yc)
    iy, ix = np.unravel_index(int(np.argmin(fe)), fe.shape)
    out = {
        "unit": "kJ/mol",
        "frames": len(pc1),
        "pc1_range": [round(float(pc1.min()), 3), round(float(pc1.max()), 3)],
        "pc2_range": [round(float(pc2.min()), 3), round(float(pc2.max()), 3)],
        "min_basin_pc1": round(float(X[iy, ix]), 3),
        "min_basin_pc2": round(float(Y[iy, ix]), 3),
        "min_G": round(float(fe[iy, ix]), 3),
    }
    if label is not None:
        try:
            lab, nlab = label(fe < FEL_BASIN_CUTOFF)
            sizes = np.bincount(lab.ravel())
            sizes = sizes[1:] if len(sizes) > 1 else np.array([])
            largest = sorted(sizes.tolist(), reverse=True)[:3]
            out["n_basins_below_2kJ"] = int(nlab)
            out["largest_basin_bin_counts"] = [int(s) for s in largest]
            total_bins = int(valid.sum())
            if largest and total_bins:
                out["largest_basin_frame_fraction"] = round(largest[0] / total_bins, 4)
        except Exception:
            pass
    return out


def _qc_block(md_dir):
    out = {}
    rows = _read_xvg(os.path.join(md_dir, "em_potential.xvg"))
    if rows:
        y = np.array([r[1] for r in rows], dtype=float)
        out["em_potential"] = {
            "unit": "kJ/mol",
            "n_points": len(y),
            "first": round(float(y[0]), 1),
            "last": round(float(y[-1]), 1),
        }
    units = {"nvt_temperature": "K", "npt_temperature": "K", "npt_pressure": "bar", "npt_density": "kg/m^3"}
    for key, fname in (("nvt_temperature", "nvt_temp.xvg"), ("npt_temperature", "npt_temp.xvg"),
                       ("npt_pressure", "npt_pres.xvg"), ("npt_density", "npt_dens.xvg")):
        rows = _read_xvg(os.path.join(md_dir, fname))
        if rows:
            y = np.array([r[1] for r in rows], dtype=float)
            out[key] = {
                "unit": units[key],
                "mean": round(float(np.mean(y)), 3),
                "std": round(float(np.std(y)), 3),
                "min": round(float(np.min(y)), 3),
                "max": round(float(np.max(y)), 3),
            }
    return out or None


def _outputs_block(analysis_dir):
    present, missing = [], []
    for level in EXPECTED_OUTPUTS:
        for name in EXPECTED_OUTPUTS[level]:
            (present if os.path.isfile(os.path.join(analysis_dir, name)) else missing).append(name)
    return {"present": present, "missing": missing}


TREND_EN = {
    "plateau": "Plateau (|slope| < 0.05 A/ns)",
    "slow_drift": "Slow drift (0.05 <= |slope| < 0.2 A/ns)",
    "evident_drift": "Evident drift (|slope| >= 0.2 A/ns)",
}


def _trend_short(cls):
    return {"plateau": "Plateau", "slow_drift": "Slow drift",
            "evident_drift": "Evident drift"}[cls]


def _esc(s):
    return _html.escape(str(s))


def _display_path(p):
    p = str(p)
    try:
        rel = os.path.relpath(p, PROJECT_DIR)
    except ValueError:
        rel = None
    if rel and not rel.startswith(".."):
        return rel.replace("\\", "/")
    return os.path.basename(p.rstrip("/\\"))


def _xvg_xunit(path):
    try:
        with open(path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if not s.startswith(("#", "@")):
                    break
                if s.startswith("@") and "xaxis" in s and "label" in s:
                    if "(ps)" in s:
                        return "ps"
                    if "(ns)" in s:
                        return "ns"
    except OSError:
        pass
    return None


def _svg_chart(path, color, x_label, y_label, scale=1.0, max_points=240):
    rows = _read_xvg(path)
    if not rows or len(rows) < 2:
        return ""
    xs = [r[0] for r in rows]
    ys = [r[1] * scale for r in rows]
    if _xvg_xunit(path) == "ps":
        xs = [v / 1000.0 for v in xs]
        if x_label.startswith("Time"):
            x_label = "Time (ns)"
    n = len(xs)
    if n > max_points:
        idx = [int(i * (n - 1) / (max_points - 1)) for i in range(max_points)]
        xs = [xs[i] for i in idx]
        ys = [ys[i] for i in idx]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    if x1 == x0:
        x1 = x0 + 1.0
    if y1 == y0:
        y1 = y0 + 1.0
    W, H, L, R, T, B = 520, 210, 62, 16, 16, 40
    pw = W - L - R
    ph = H - T - B

    def px(v):
        return L + (v - x0) / (x1 - x0) * pw

    def py(v):
        return T + (1.0 - (v - y0) / (y1 - y0)) * ph

    pts = " ".join(f"{px(a):.1f},{py(b):.1f}" for a, b in zip(xs, ys))

    def xfmt(v):
        return f"{v:.0f}" if abs(x1 - x0) >= 10 else f"{v:.2f}"

    def yfmt(v):
        return f"{v:.1f}" if abs(y1 - y0) >= 1 else f"{v:.3f}"

    parts = []
    parts.append(f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
                 f'xmlns="http://www.w3.org/2000/svg" role="img">')
    parts.append(f'<rect x="{L}" y="{T}" width="{pw}" height="{ph}" fill="#fdfdfe" stroke="#e5e7eb"/>')
    for frac in (0.25, 0.5, 0.75):
        gy = T + ph * frac
        parts.append(f'<line x1="{L}" y1="{gy:.1f}" x2="{L + pw}" y2="{gy:.1f}" '
                     f'stroke="#eef1f5" stroke-width="1"/>')
    parts.append(f'<polyline points="{pts}" fill="none" stroke="{color}" '
                 f'stroke-width="1.6" stroke-linejoin="round"/>')
    parts.append(f'<text x="{L - 6}" y="{T + 4}" font-size="10.5" fill="#6b7280" '
                 f'text-anchor="end">{yfmt(y1)}</text>')
    parts.append(f'<text x="{L - 6}" y="{T + ph / 2 + 3}" font-size="10.5" fill="#6b7280" '
                 f'text-anchor="end">{yfmt((y0 + y1) / 2)}</text>')
    parts.append(f'<text x="{L - 6}" y="{T + ph + 3}" font-size="10.5" fill="#6b7280" '
                 f'text-anchor="end">{yfmt(y0)}</text>')
    parts.append(f'<text x="{L}" y="{H - 10}" font-size="10.5" fill="#6b7280" '
                 f'text-anchor="middle">{xfmt(x0)}</text>')
    parts.append(f'<text x="{L + pw / 2}" y="{H - 10}" font-size="10.5" fill="#6b7280" '
                 f'text-anchor="middle">{xfmt((x0 + x1) / 2)}</text>')
    parts.append(f'<text x="{L + pw}" y="{H - 10}" font-size="10.5" fill="#6b7280" '
                 f'text-anchor="middle">{xfmt(x1)}</text>')
    parts.append(f'<text x="{W / 2}" y="{H - 26}" font-size="11" fill="#374151" '
                 f'text-anchor="middle">{_esc(x_label)}</text>')
    parts.append(f'<text x="14" y="{T + ph / 2}" font-size="11" fill="#374151" '
                 f'text-anchor="middle" transform="rotate(-90 14 {T + ph / 2})">{_esc(y_label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def _card(label, value, sub=""):
    sub_html = f'<span class="sub">{_esc(sub)}</span>' if sub else ""
    return (f'<div class="card"><div class="card-label">{_esc(label)}</div>'
            f'<div class="card-value">{_esc(value)}{sub_html}</div></div>')


def _video_tag():
    if not os.path.isfile(LOGO_VIDEO):
        return ""
    with open(LOGO_VIDEO, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("ascii")
    return (f'<video class="logo" autoplay muted loop playsinline '
            f'aria-hidden="true" title="DockMD logo" '
            f'src="data:video/mp4;base64,{b64}"></video>')


FEL_LINKS = {
    "npt": [
        ("04_FEL_Plot_Advanced.png",
         "04_FEL_Plot_Advanced.png — Gaussian-smoothed 2D free energy landscape (NPT)"),
        ("06_FEL_3D_Interactive.html",
         "06_FEL_3D_Interactive.html — interactive 3D free energy landscape (NPT)"),
    ],
    "md": [
        ("09_md_FEL_Plot_Advanced.png",
         "09_md_FEL_Plot_Advanced.png — Gaussian-smoothed 2D free energy landscape (production)"),
        ("11_md_FEL_3D_Interactive.html",
         "11_md_FEL_3D_Interactive.html — interactive 3D free energy landscape (production)"),
    ],
}


def _fel_links_html(analysis_dir, level):
    items = []
    for fname, caption in FEL_LINKS.get(level, []):
        if not os.path.isfile(os.path.join(analysis_dir, fname)):
            continue
        hint = ("open in a browser; drag to rotate, scroll to zoom"
                if fname.endswith(".html") else "static image")
        items.append(f'<li><a href="{_esc(fname)}" target="_blank">{_esc(caption)}</a> '
                     f'<span class="hint">— {hint}</span></li>')
    if not items:
        return ""
    title = "NPT stage" if level == "npt" else "production stage"
    return ('<h3 style="font-size:14px;margin:18px 0 8px">Free energy landscape '
            f'files ({title})</h3><ul class="file-links">' + "".join(items) + "</ul>")


def _level_section_html(level, blk, charts, analysis_dir=""):
    if not any(blk.values()):
        return ""
    title = ("NPT Equilibration Trajectory Dynamics (npt level)"
             if level == "npt" else "Production Trajectory Dynamics (md.xtc)")
    h = [f"<h2>{_esc(title)}</h2>"]
    cards = []
    r = blk.get("rmsd")
    if r:
        cards.append(_card("Backbone RMSD, final window",
                           f"{r['late_mean']} ± {r['late_std']} Å",
                           f"initial {r['initial_mean']} Å"))
        cards.append(_card("RMSD trend class", _trend_short(r["trend_class"]),
                           f"max {r['max']} Å @ {r['max_time_ns']} ns"))
        cards.append(_card("Trajectory length", f"{r['time_span_ns']} ns",
                           f"{r['frames']} frames"))
    g = blk.get("gyrate")
    if g:
        cards.append(_card("Radius of gyration Rg, final window",
                           f"{g['late_mean']} nm",
                           f"change {g['change_percent']}%"))
    f_ = blk.get("rmsf")
    if f_:
        cards.append(_card("Residue RMSF, mean / max",
                           f"{f_['mean']} / {f_['max']} Å",
                           f"{f_['n_above_mean_plus_2sd']} residues above mean + 2σ"))
    hbn = blk.get("hbnum")
    if hbn:
        cards.append(_card("Protein-ligand H-bonds",
                           f"mean {hbn['mean']} per frame",
                           f"max {hbn['max']}; frames with ≥ 1 H-bond: "
                           f"{hbn['fraction_frames_with_at_least_one'] * 100:.1f}%"))
    fe = blk.get("fel")
    if fe:
        cards.append(_card("FEL global minimum", f"G ≈ {fe['min_G']} kJ/mol",
                           f"(PC1 ≈ {fe['min_basin_pc1']}, PC2 ≈ {fe['min_basin_pc2']})"))
    h.append('<div class="cards">' + "".join(cards) + "</div>")
    if f_ and f_.get("top_flexible_residues"):
        top = ", ".join(f"res {t['residue']} ({t['rmsf_A']} Å)"
                        for t in f_["top_flexible_residues"])
        h.append(f'<ul class="plain"><li>Most flexible residues: {_esc(top)}. '
                 f'<span style="color:#6b7280;font-size:12px">(Residue numbers are '
                 f'GROMACS selection indices and may be offset from the PDB numbering.)'
                 f'</span></li></ul>')
    figs = []
    names = [("rmsd", "Protein Backbone RMSD"),
             ("gyrate", "Radius of Gyration (Rg)"),
             ("rmsf", "Residue RMSF"),
             ("hbnum", "Protein-Ligand H-bonds")]
    for key, caption in names:
        svg = charts.get(key, "")
        if svg:
            figs.append(f"<figure>{svg}<figcaption>{_esc(caption)}</figcaption></figure>")
    if figs:
        h.append('<div class="charts">' + "".join(figs) + "</div>")
    if fe and "n_basins_below_2kJ" in fe:
        h.append(f'<ul class="plain"><li>Distinct basins below '
                 f'{FEL_BASIN_CUTOFF} kJ/mol: {fe["n_basins_below_2kJ"]} '
                 f'(largest basin occupies '
                 f'{fe.get("largest_basin_frame_fraction", "NA")} of the populated bins)</li></ul>')
    links = _fel_links_html(analysis_dir, level)
    if links:
        h.append(links)
    return "".join(h)


def _qc_section_html(qc):
    if not qc:
        return ""
    h = ["<h2>Equilibration QC (NPT stage)</h2>"]
    cards = []
    em = qc.get("em_potential")
    if em:
        cards.append(_card("Energy minimisation, potential",
                           f"reduced to {em['last']} kJ/mol",
                           f"initial {em['first']}, {em['n_points']} recorded points"))
    for key, label in (("nvt_temperature", "NVT temperature"),
                       ("npt_temperature", "NPT temperature"),
                       ("npt_pressure", "NPT pressure"),
                       ("npt_density", "NPT density")):
        blk = qc.get(key)
        if blk:
            cards.append(_card(label, f"{blk['mean']} ± {blk['std']} {blk['unit']}",
                               f"range {blk['min']} – {blk['max']}"))
    h.append('<div class="cards">' + "".join(cards) + "</div>")
    h.append('<ul class="plain"><li>Mean temperature close to the target and mean '
             'density close to the solvent reference, with limited fluctuations, '
             'usually indicate adequate equilibration; large instantaneous NPT '
             'pressure fluctuations are expected and not used as a criterion.</li></ul>')
    return "".join(h)


def _outputs_section_html(out):
    if not out:
        return ""
    h = ['<h2>Output files</h2><table><tr><th>Status</th><th>File</th></tr>']
    for name in out.get("present", []):
        h.append(f'<tr><td class="ok">Generated</td><td>{_esc(name)}</td></tr>')
    for name in out.get("missing", []):
        h.append(f'<tr><td class="miss">Missing</td><td>{_esc(name)}'
                 f'<span style="color:#6b7280"> (corresponding analysis level '
                 f'has not been run)</span></td></tr>')
    h.append("</table>")
    return "".join(h)


def _rmsd_interpretation(r):
    if not r:
        return ""
    parts = [f"Backbone RMSD averaged {r['initial_mean']} Å over the initial 5% of "
             f"the trajectory and {r['late_mean']} ± {r['late_std']} Å over the final "
             f"20% ({r['late_window_ns'][0]}–{r['late_window_ns'][1]} ns); the overall "
             f"maximum was {r['max']} Å at {r['max_time_ns']} ns."]
    slope = r["tail_slope_A_per_ns"]
    tc = r["trend_class"]
    if r.get("short_track"):
        return (f"The trajectory spans only {r['time_span_ns']} ns and serves as an "
                f"equilibration spot-check: a trend classification of "
                f"{_trend_short(tc).lower()} (slope {slope} Å/ns) is not reliable at "
                f"this length and only confirms that the system evolves normally; "
                f"stability conclusions should be based on the production (md) level. "
                f"Final-window mean: {r['late_mean']} ± {r['late_std']} Å.")
    if tc == "plateau":
        parts.append(f"The tail slope of {slope} Å/ns (|slope| < {RMSD_PLATEAU}) "
                     f"indicates a plateau; the deviation from the initial "
                     f"conformation is essentially stable.")
    elif tc == "slow_drift":
        parts.append(f"The tail slope of {slope} Å/ns falls in the slow-drift band "
                     f"({RMSD_PLATEAU}–{RMSD_DRIFT} Å/ns): the system is still "
                     f"relaxing gradually, and windowed means rather than "
                     f"single-point values are recommended when reporting statistics.")
    else:
        parts.append(f"The tail slope of {slope} Å/ns (>= {RMSD_DRIFT} Å/ns) "
                     f"indicates evident drift: the trajectory has not reached a "
                     f"plateau, and stability conclusions at this trajectory length "
                     f"should be considered preliminary.")
    return "".join(parts)


def _rg_interpretation(g):
    if not g:
        return ""
    chg = g["change_percent"]
    if chg < -1.0:
        direction = "an overall compaction"
    elif chg > 1.0:
        direction = "an overall expansion"
    else:
        direction = "an essentially unchanged overall size"
    return (f"The radius of gyration changed from {g['early_mean']} nm (early "
            f"window) to {g['late_mean']} ± {g['late_std']} nm (final window; "
            f"{chg}%), i.e., {direction}. A stable Rg once the backbone RMSD "
            f"plateaus usually indicates no global collapse or swelling of the "
            f"fold; localised changes should be localised via the RMSF profile.")


def _rmsf_interpretation(f_):
    if not f_:
        return ""
    top = ", ".join(f"res {t['residue']} ({t['rmsf_A']} Å)"
                    for t in f_["top_flexible_residues"][:5])
    s = (f"Per-residue RMSF averaged {f_['mean']} ± {f_['std']} Å "
         f"(maximum {f_['max']} Å). The most flexible residues are: {top}. "
         f"{f_['n_above_mean_plus_2sd']} residues exceed the mean + 2σ.")
    if f_["max"] > 3.0:
        s += (" The markedly high maximum typically arises from the N/C termini or "
              "flexible loop regions; if it corresponds to a loop rebuilt during "
              "structure preparation, the flexibility reflects the modelled "
              "geometry rather than an experimentally observed conformation, and "
              "the two cases should be distinguished.")
    s += (" (Residue numbers are GROMACS selection indices and may be offset from "
          "the source PDB numbering.)")
    return s


def _hb_interpretation(h):
    if not h:
        return ""
    occ = h["fraction_frames_with_at_least_one"] * 100
    s = (f"The protein-ligand hydrogen-bond count averaged {h['mean']} per frame "
         f"(maximum {h['max']}); frames containing at least one H-bond account for "
         f"{occ:.1f}% (0.35 nm heavy-atom donor-acceptor criterion).")
    if occ >= 50.0:
        s += (" H-bond contacts persist in the majority of frames, indicating a "
              "fairly stable direct H-bond contribution at the interface.")
    elif occ >= 10.0:
        s += (" H-bonds occur intermittently; binding is more likely dominated by "
              "hydrophobic and van der Waals contacts, with hydrogen bonds "
              "contributing in a transient, on-off fashion (descriptive statement, "
              "not a causal conclusion).")
    else:
        s += (" The H-bond occupancy is low; interfacial contacts are dominated by "
              "non-polar interactions (descriptive statement, not a causal "
              "conclusion).")
    return s


def _fel_interpretation(fe):
    if not fe:
        return ""
    if fe["frames"] < FEL_MIN_FRAMES:
        return (f"The free energy landscape was projected from {fe['frames']} "
                f"frames, below the {FEL_MIN_FRAMES}-frame threshold; basin counts "
                f"and positions are not quantitatively meaningful at this sampling "
                f"depth, and only the approximate location of the global minimum "
                f"should be referenced.")
    s = (f"The global minimum of the free energy surface lies at "
         f"G ≈ {fe['min_G']} kJ/mol, at (PC1 ≈ {fe['min_basin_pc1']}, "
         f"PC2 ≈ {fe['min_basin_pc2']}).")
    n = fe.get("n_basins_below_2kJ")
    if n is not None:
        frac = fe.get("largest_basin_frame_fraction")
        s += (f" {n} distinct basins lie below {FEL_BASIN_CUTOFF} kJ/mol"
              + (f", the largest accounting for {frac} of the populated bins."
                 if frac is not None else "."))
        if n <= 2:
            s += (" The small number of basins suggests conformational sampling "
                  "concentrated in a few substates.")
        else:
            s += (" The larger number of basins suggests multiple mutually "
                  "accessible substates; interconversion frequencies at this "
                  "trajectory length should not be over-interpreted.")
    return s


def _qc_interpretation(qc, targets):
    if not qc:
        return ""
    parts = []
    em = qc.get("em_potential")
    if em:
        parts.append(f"<li>Energy minimisation: the potential decreased from "
                     f"{em['first']} to {em['last']} kJ/mol over "
                     f"{em['n_points']} recorded points, indicating normal "
                     f"convergence.</li>")
    t_target = targets.get("temperature", 300.0)
    p_target = targets.get("pressure", 1.0)
    src = targets.get("source", "defaults")
    checks = (("nvt_temperature", "NVT temperature", t_target, 5.0),
              ("npt_temperature", "NPT temperature", t_target, 5.0),
              ("npt_density", "NPT density", 1000.0, 30.0))
    for key, label, target, tol in checks:
        blk = qc.get(key)
        if blk:
            ok = abs(blk["mean"] - target) <= tol
            verdict = ("within tolerance" if ok else
                       "deviates from the target; check temperature coupling or "
                       "parameter settings")
            parts.append(f"<li>{label}: mean {blk['mean']} ± {blk['std']} "
                         f"{blk['unit']} (target {target:.0f} ± {tol:.0f}, range "
                         f"{blk['min']} – {blk['max']}): <b>{verdict}</b>.</li>")
    if targets.get("pressure") is not None:
        parts.append(f"<li>Target temperature and pressure were taken from {src} "
                     f"(ref_t = {t_target:g} K, ref_p = {p_target:g} bar); the NPT "
                     f"density target of 1000 kg/m³ assumes an aqueous solvent and "
                     f"must be revised manually for non-aqueous systems.</li>")
    pres = qc.get("npt_pressure")
    if pres:
        parts.append(f"<li>NPT pressure: mean {pres['mean']} ± {pres['std']} bar "
                     f"(range {pres['min']} – {pres['max']}). Large instantaneous "
                     f"pressure fluctuations (hundreds of bar) are normal for "
                     f"small NPT ensembles and are not used as a criterion.</li>")
    return "<ul class=\"plain\">" + "".join(parts) + "</ul>"


def _interpretation_section_html(data):
    blocks = []
    qc_txt = _qc_interpretation(data.get("qc"), data.get("mdp_targets") or {})
    if qc_txt:
        blocks.append(("Equilibration QC", qc_txt))
    for lv, label in (("npt", "NPT equilibration trajectory (npt level)"),
                      ("md", "Production trajectory (md level)")):
        blk = data["levels"].get(lv) or {}
        items = []
        for key, fn in (("rmsd", _rmsd_interpretation), ("gyrate", _rg_interpretation),
                        ("rmsf", _rmsf_interpretation), ("hbnum", _hb_interpretation),
                        ("fel", _fel_interpretation)):
            txt = fn(blk.get(key))
            if txt:
                items.append(f"<li>{txt}</li>")
        if items:
            blocks.append((label, "<ul class=\"plain\">" + "".join(items) + "</ul>"))
    if not blocks:
        return ""
    h = ['<h2>Interpretation (auto-generated)</h2>',
         '<div class="notice">The interpretations below are computed at report '
         'time by the deterministic statistics script directly from the xvg data; '
         'they share the same source as the charts, are reproducible, and require '
         'no manual backfill. They describe statistically observable features '
         'only and do not constitute causal inference.</div>']
    for title, body in blocks:
        h.append(f'<div class="interp"><h3>{_esc(title)}</h3>'
                 f'<div class="interp-body">{body}</div></div>')
    return "".join(h)


def _ai_section_html(data):
    meta = data.get("ai_interpretation")
    if not meta or meta.get("status") not in ("ok", "degraded"):
        return ""
    dropped = meta.get("dropped_sentences", 0)
    warn = ""
    if meta["status"] == "degraded":
        warn = (f'<div class="notice-amber">{dropped} sentences failed the '
                f'numeric verification and were discarded automatically; the '
                f'remaining text is consistent with the JSON statistics.</div>')
    h = [f'<h2>AI-Assisted Interpretation ({_esc(meta["model"])})</h2>',
         '<div class="notice-amber">This section was generated by a large '
         'language model from 12_AI_Report_Data.json (temperature = 0). Every '
         'numeric statement was programmatically verified against the JSON '
         'statistics; the narrative only describes statistically observable '
         'features. Formal conclusions require manual review.</div>',
         warn,
         f'<div class="interp"><div class="interp-body">{meta["text_html"]}</div></div>']
    return "".join(h)


def _report_html(data, charts):
    levels_used = [lv for lv in ("npt", "md") if any(data["levels"].get(lv, {}).values())]
    h = ['<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
         f"<title>DockMD Automated Analysis Report - {_esc(data['system'])}</title>",
         f"<style>{CSS}</style></head><body><div class=\"wrap\">",
         '<header><div class="title"><h1>DockMD Automated Analysis Report</h1>'
         f"<div class=\"meta\">System: <b>{_esc(data['system'])}</b><br>"
         f"Generated: {_esc(data['generated_at'])}<br>"
         f"MD intermediate data: {_esc(_display_path(data['md_dir']))}<br>"
         f"Output directory: {_esc(_display_path(data['analysis_dir']))}<br>"
         f"Levels covered: {_esc(', '.join(levels_used) if levels_used else 'no data available')}</div></div>",
         _video_tag(),
         "</header>"]
    h.append('<div class="notice">This report is generated deterministically '
             'from the xvg intermediate files by the reporting script; it '
             'provides descriptive summaries only and does not imply causal '
             'relationships. The AI-Assisted Interpretation section, when '
             'present, is for reference only and requires manual review before '
             'use in formal materials.</div>')
    h.append(_qc_section_html(data.get("qc")))
    for lv in ("npt", "md"):
        h.append(_level_section_html(lv, data["levels"].get(lv, {}) or {},
                                     charts.get(lv, {}), data["analysis_dir"]))
    h.append(_outputs_section_html(data.get("outputs")))
    h.append(_interpretation_section_html(data))
    h.append(_ai_section_html(data))
    h.append("</div></body></html>")
    return "".join(h)


def _level_charts(md_dir, level):
    charts = {}
    for key, fname, color, xl, yl, scale in (
            ("rmsd", f"{level}_rmsd.xvg", "#1f77b4", "Time (ns)", "RMSD (Å)", 10.0),
            ("gyrate", f"{level}_gyrate.xvg", "#2ca02c", "Time (ns)", "Rg (nm)", 1.0),
            ("rmsf", f"{level}_rmsf.xvg", "#d62728", "Residue", "RMSF (Å)", 10.0),
            ("hbnum", f"{level}_hbnum.xvg", "#9467bd", "Time (ns)", "H-bonds", 1.0)):
        svg = _svg_chart(os.path.join(md_dir, fname), color, xl, yl, scale)
        if svg:
            charts[key] = svg
    return charts


def build_report(md_dir, levels=("npt", "md"), analysis_dir=None, use_ai=False):
    if analysis_dir is None:
        analysis_dir = os.path.join(ANALYSIS_ROOT, os.path.basename(md_dir))
    os.makedirs(analysis_dir, exist_ok=True)
    data = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "system": os.path.basename(md_dir),
        "md_dir": md_dir,
        "analysis_dir": analysis_dir,
        "levels": {},
        "qc": _qc_block(md_dir) if "npt" in levels else None,
        "mdp_targets": _mdp_targets(md_dir) if "npt" in levels else None,
        "outputs": _outputs_block(analysis_dir),
    }
    charts = {}
    for level in levels:
        data["levels"][level] = {
            "rmsd": _rmsd_block(md_dir, level),
            "gyrate": _gyrate_block(md_dir, level),
            "rmsf": _rmsf_block(md_dir, level),
            "hbnum": _hbnum_block(md_dir, level),
            "fel": _fel_block(md_dir, level),
        }
        charts[level] = _level_charts(md_dir, level)
    if use_ai:
        from ai import interpret
        try:
            data["ai_interpretation"] = interpret.generate_interpretation(data, cfg)
        except Exception as e:
            data["ai_interpretation"] = {"model": cfg.LLM_MODEL, "status": "failed",
                                         "dropped_sentences": 0, "text_html": "",
                                         "error": f"{type(e).__name__}: {e}"}
        status = data["ai_interpretation"]["status"]
        if status == "ok":
            print("[OK] AI interpretation generated (all numbers verified).")
        elif status == "degraded":
            print(f"[OK] AI interpretation generated "
                  f"({data['ai_interpretation']['dropped_sentences']} sentences "
                  f"dropped by numeric verification).")
        else:
            print(f"[WARN] AI interpretation unavailable: "
                  f"{data['ai_interpretation']['error']}")
            print("       Report continues with deterministic content only.")
    json_path = os.path.join(analysis_dir, "12_AI_Report_Data.json")
    html_path = os.path.join(analysis_dir, "13_AI_Report.html")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(_report_html(data, charts))
    return {"report_data": json_path, "report_html": html_path, "data": data}


def main(level="all", use_ai=False):
    state = load_state()
    md_dir = state.get("md_dir")
    if not md_dir or not os.path.isdir(md_dir):
        sys.exit("!! MD working directory not found in state.json. Run the MD steps first.")
    levels = ("npt", "md") if level == "all" else (level,)
    res = build_report(md_dir, levels, use_ai=use_ai)
    print(f"\n[OK] Report data: {res['report_data']}")
    print(f"[OK] Report HTML: {res['report_html']}")
    print("Next: open 13_AI_Report.html in a browser. Interpretation text is "
          "generated deterministically from the xvg data at report time.")
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", default="all", choices=["all", "npt", "md"])
    ap.add_argument("--ai", action="store_true",
                    help="append LLM interpretation section (needs API key env)")
    ns = ap.parse_args()
    main(level=ns.level, use_ai=ns.ai)
