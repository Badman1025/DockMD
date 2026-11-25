import json
import os
import re
import urllib.error
import urllib.request

NUM_RE = re.compile(r"-?\d+\.\d+|-?\d+")

SYSTEM_PROMPT = (
    "You are a molecular dynamics analysis assistant writing the AI-assisted "
    "interpretation section of a formal MD analysis report. You receive a JSON "
    "summary of deterministic statistics computed from GROMACS xvg data "
    "(backbone RMSD, radius of gyration Rg, per-residue RMSF, protein-ligand "
    "hydrogen bonds, free energy landscape FEL, equilibration QC). Write in "
    "English using precise, professional MD terminology and a formal, "
    "impersonal academic register. Strictly observe these rules:\n"
    "1. Use only numbers that appear in the JSON. Never invent, estimate, or "
    "introduce any number from outside the JSON; every numeric statement must "
    "match the JSON exactly, and every number must carry its unit "
    "(A, nm, ns, kJ/mol, K, bar, kg/m^3).\n"
    "2. Descriptive interpretation and cross-metric consistency review only. "
    "No causal inference. Do not judge ligand binding affinity, residence, or "
    "therapeutic efficacy. Phrase uncertain statements with appropriate "
    "hedging ('suggests', 'is consistent with', 'indicates') rather than "
    "definitive claims.\n"
    "3. Style: concise declarative sentences; no first person, no "
    "colloquialisms, no exclamations, no rhetorical questions; use standard "
    "terms (plateau, convergence, drift, occupancy, sampling depth, substate, "
    "basin, equilibration).\n"
    "4. Structure the output into four subsections in this exact order, each "
    "starting on its own line with exactly one of these labels followed by a "
    "colon: \"Equilibration stage (NPT):\", \"Production stage:\", "
    "\"Cross-stage comparison:\", \"Summary and caveats:\". Under each label "
    "use concise bullet points (one sentence each) or a short paragraph.\n"
    "5. \"Equilibration stage (NPT)\" must describe only the statistics under "
    "levels.npt (RMSD, Rg, RMSF, hydrogen bonds, FEL).\n"
    "6. \"Production stage\" must describe the statistics under levels.md if "
    "that key is present and non-empty. If no production-level data is "
    "present in the JSON, state that explicitly and recommend running the "
    "production MD step next; do not invent production-level statements.\n"
    "7. \"Cross-stage comparison\" must explicitly compare the short NPT "
    "trajectory with the long production trajectory: RMSD magnitude and trend "
    "class, Rg change, FEL basin structure, and hydrogen-bond occupancy. If "
    "no production-level data exists, state that a cross-stage comparison is "
    "not yet possible.\n"
    "8. \"Summary and caveats\" must close with a brief overall summary plus "
    "data caveats (sampling depth, residue-numbering offset, QC tolerance).\n"
    "9. When referring to an RMSD trend classification, quote it verbatim in "
    "its canonical form (\"Plateau\", \"Slow drift\", \"Evident drift\") as "
    "given by the trend_class field. Never weaken or strengthen a "
    "deterministic classification (e.g. do not call an \"Evident drift\" a "
    "\"minor drift\").\n"
    "10. A drift observed in the short NPT trajectory reflects ongoing "
    "structural relaxation during equilibration and is expected behavior; do "
    "not describe it as an equilibration failure or an equilibration issue "
    "unless a QC metric is outside tolerance.\n"
    "11. Do not speculate beyond the JSON (e.g., residue function or disease "
    "associations); where structural context is required, state that it "
    "'requires manual inspection'. Residue numbers are GROMACS selection "
    "indices and may be offset from the PDB numbering.\n"
    "12. Keep the total length under 350 words."
)

SMALL_INT_MAX = 20
REL_TOL = 0.05
ABS_TOL = 0.011


def _collect_values(obj, out):
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        out.append(float(obj))
    elif isinstance(obj, dict):
        for v in obj.values():
            _collect_values(v, out)
    elif isinstance(obj, list):
        for v in obj:
            _collect_values(v, out)


def _number_ok(num, values):
    for v in values:
        if abs(num - v) <= max(REL_TOL * abs(v), ABS_TOL):
            return True
    return False


def _split_sentences(text):
    sentences = []
    last = 0
    for m in re.finditer(r"[。！？!?；;\n]|\.(?=\s)", text):
        seg = text[last:m.end()]
        sentences.append(seg)
        last = m.end()
    tail = text[last:]
    sentences.append(tail)
    return sentences


def verify_numbers(text, data):
    values = []
    _collect_values(data, values)
    sentences = _split_sentences(text)
    kept, dropped = [], 0
    for s in sentences:
        if not s.strip():
            kept.append(s)
            continue
        nums = [float(n) for n in NUM_RE.findall(s)]
        bad = [n for n in nums
               if not (n.is_integer() and abs(n) <= SMALL_INT_MAX)
               and not _number_ok(n, values)]
        if bad:
            dropped += 1
            continue
        kept.append(s)
    return "".join(kept), dropped


SECTION_LABELS = ("Equilibration stage (NPT)", "Production stage",
                  "Cross-stage comparison", "Summary and caveats")

_LABEL_RE = re.compile(r"^(" + "|".join(re.escape(l) for l in SECTION_LABELS)
                       + r")\s*[:：]\s*(.*)$")


def _markdown_to_html(text):
    import html as _html
    lines = []
    for raw in text.splitlines():
        s = raw.strip()
        if not s:
            continue
        s = _html.escape(s)
        m = _LABEL_RE.match(s)
        if m:
            lines.append("<p style='margin:10px 0 4px;font-weight:600;"
                         "color:#1d4ed8'>" + m.group(1) + "</p>")
            if m.group(2).strip():
                lines.append(f"<p style='margin:4px 0'>{m.group(2).strip()}</p>")
        elif s.startswith("- ") or s.startswith("* "):
            lines.append(f"<li>{s[2:]}</li>")
        elif re.match(r"^\d+\. ", s):
            stripped = re.sub(r"^\d+\. ", "", s)
            lines.append(f"<li>{stripped}</li>")
        else:
            lines.append(f"<p style='margin:4px 0'>{s}</p>")
    out = []
    in_list = False
    for item in lines:
        if item.startswith("<li>"):
            if not in_list:
                out.append("<ul class=\"plain\">")
                in_list = True
            out.append(item)
        else:
            if in_list:
                out.append("</ul>")
                in_list = False
            out.append(item)
    if in_list:
        out.append("</ul>")
    return "".join(out)


def generate_interpretation(data, cfg):
    env_name = (cfg.LLM_API_KEY_ENV or "").strip()
    api_key = os.environ.get(env_name, "")
    meta = {"model": cfg.LLM_MODEL, "status": "failed", "dropped_sentences": 0,
            "text_html": "", "error": "", "key_source": f"env {env_name}"}
    if not api_key:
        meta["error"] = (f"missing API key: run 'export {env_name}=sk-...' "
                         f"before generating the report")
        return meta
    payload = {
        "model": cfg.LLM_MODEL,
        "temperature": 0,
        "messages": [
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user",
             "content": "Below is the DockMD deterministic statistics summary "
                        "JSON. Interpret it according to the system rules:\n"
                        + json.dumps(data, ensure_ascii=False)},
        ],
    }
    url = cfg.LLM_BASE_URL.rstrip("/") + "/chat/completions"
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json",
                 "Authorization": f"Bearer {api_key}"},
        method="POST")
    try:
        with urllib.request.urlopen(req, timeout=180) as resp:
            body = json.loads(resp.read().decode("utf-8"))
        text = body["choices"][0]["message"]["content"].strip()
    except (urllib.error.URLError, urllib.error.HTTPError, KeyError,
            IndexError, ValueError, OSError) as e:
        meta["error"] = f"{type(e).__name__}: {e}"
        return meta
    kept, dropped = verify_numbers(text, data)
    meta["dropped_sentences"] = dropped
    if not kept.strip():
        meta["status"] = "rejected"
        meta["error"] = "all sentences failed numeric verification"
        return meta
    meta["status"] = "ok" if dropped == 0 else "degraded"
    meta["text_html"] = _markdown_to_html(kept)
    return meta
