"""Build REPORT.md from report_template.md, the result files and the source code.
usage: python build_report.py        (run inside detected_event_extrapolation/)"""
import re, json, os, numpy as np
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
P = lambda *a: os.path.join(HERE, *a)

# ---------------- Clifford table, parsed from results.txt
def clifford_table():
    blocks, cur = {}, None
    for line in open(P("results.txt")):
        m = re.match(r"== (\w): (.*?):? n=\d+ k=(\d+) reps=\d+ p=([\d.]+)", line)
        if m: cur = m.group(1); blocks[cur] = dict(title=m.group(2), k=m.group(3), p=m.group(4), med={}, rms={}); continue
        if cur is None: continue
        m = re.match(r"\s+(\w+)\s+median=\s*([\d.eE+-]+)", line)
        if m: blocks[cur]["med"][m.group(1)] = float(m.group(2))
        m = re.match(r"\s+bias\[(.*?)\s*\] mean=([\d.eE+-]+)\s+rms=([\d.eE+-]+)", line)
        if m: blocks[cur]["rms"][m.group(1)] = float(m.group(3))
        if line.startswith("== scaling"): cur = None
    names = {"A": "A: checks wrap whole circuit", "B": "B: same, 3× noise", "C": "C: only 2 checks", "D": "D: checks on 2nd half", "E": "E: checks on 1st half"}
    keys = [k for k in "ABCDE" if k in blocks]
    out = ["| | " + " | ".join(names[k] for k in keys) + " |", "|---|" + "---|" * len(keys)]
    row = lambda lab, f: out.append(f"| {lab} | " + " | ".join(f(blocks[k]) for k in keys) + " |")
    row("checks / error rate p", lambda b: f"{b['k']} / {b['p']}")
    row("total fault rate Λ", lambda b: f"{b['med']['Lam']:.2f}")
    row("coverage / acceptance", lambda b: f"{b['med']['cov']*100:.0f}% / {b['med']['alpha']*100:.0f}%")
    for key, lab in [("raw", "unmitigated"), ("ED", "error detection only"), ("ED+ext(generic r)", "extrapolated, coverage ratio"),
                     ("ED+ext(O-specific r)", "extrapolated, observable-specific ratio"), ("ED+ext(O-spec r)+2nd", "+ pair term from rejected shots"),
                     ("ED+model rescale", "model rescaling, exact model"), ("ED+model rescale, scale off", "model rescaling, model 30% high"),
                     ("ED+ext(O-spec r), scale off", "extrapolated, model 30% high")]:
        row(lab, lambda b, key=key: f"{b['rms'][key]*100:.2f}%" if b['rms'][key] >= 5e-4 else f"{b['rms'][key]*100:.3f}%")
    row("cost: plain PEC", lambda b: f"{b['med']['G_pec']:.1f}")
    row("cost: 1st-order PEC + detection", lambda b: f"{b['med']['G_edpec']:.1f}")
    row("cost: two-point extrapolation", lambda b: f"{b['med']['G_ext_generic']:.1f}")
    return "\n".join(out) + "\n\nBias rows are rms over the 40 instances. Cost rows are $N\\,\\mathrm{Var}$ proxies for a generic $\\pm1$ observable (medians)."

# ---------------- non-Clifford tables, from the .jsonl files
recs = [json.loads(l) for fn in ("flagship.jsonl", "generic_lo.jsonl", "generic_hi.jsonl") for l in open(P("nonclifford", fn))]
G = defaultdict(list)
for r in recs: G[r["tag"]].append(r)
ROWS = [("raw", "unmitigated"), ("ED", "error detection only"), ("ext, coverage ratio", "extrapolated, coverage ratio"),
        ("ext, oracle 1st-order ratio", "extrapolated, exact 1st-order ratio"), ("ext, oracle ratio + pair term", "… + in-situ pair term"),
        ("ext, site-pooled training ratio", "**extrapolated, Clifford-trained ratio**"), ("CDR on ED (pooled slope)", "Clifford regression on ED data"),
        ("PEC 1st order (exact model)", "1st-order PEC, exact model"), ("PEC 1st order (model 30% high)", "1st-order PEC, model 30% high"),
        ("u-amp ZNE, exp fit (G=1,2,3)", "ZNE, undetected faults amplified"), ("ZNE on ED, exp fit (G=1,2,3)", "ZNE on ED, all noise amplified"),
        ("ZNE raw, exp fit (G=1,2,3)", "ZNE, no checks")]
def cols(p): return [(f"flag d=6 p={p}", "anc", "Paper angles, 6 steps"), (f"flag d=6 p={p}", "parity", "Paper angles + parity"),
                     (f"generic p={p}", "anc", "Generic ensemble"), (f"generic p={p}", "parity", "Generic + parity")]
def rms(tag, cfg, key):
    n = G[tag][0]["n"]; B = np.array([np.array(r[cfg]["est"][key])[:n] - np.array(r["ideal"])[:n] for r in G[tag]]); return np.sqrt((B ** 2).mean())
def nc_table(p):
    C = cols(p); out = ["| | " + " | ".join(c[2] for c in C) + " |", "|---|" + "---|" * len(C)]
    out.append("| total fault rate Λ | " + " | ".join(f"{np.mean([r[c]['Ld']+r[c]['Lu'] for r in G[t]]):.2f}" for t, c, _ in C) + " |")
    out.append("| coverage / acceptance | " + " | ".join(f"{np.mean([r[c]['coverage'] for r in G[t]])*100:.0f}% / {np.mean([r[c]['alpha'] for r in G[t]])*100:.0f}%" for t, c, _ in C) + " |")
    out.append("| mean ideal magnitude | " + " | ".join(f"{np.mean(np.abs([r['ideal'][:r['n']] for r in G[t]])):.2f}" for t, c, _ in C) + " |")
    for key, lab in ROWS: out.append(f"| {lab} | " + " | ".join(f"{rms(t, c, key):.4f}" for t, c, _ in C) + " |")
    return "\n".join(out)
def nc_ratios():
    out = ["| Scenario | coverage ratio | exact 1st-order (site mean) | Clifford-trained (pooled) |", "|---|---|---|---|"]
    for p in ("0.004", "0.012"):
        for t, c, lab in cols(p):
            out.append(f"| {lab}, p = {p} | {np.mean([r[c]['r_cov'] for r in G[t]]):.2f} | {np.mean([np.mean(r[c]['r1'][:r['n']]) for r in G[t]]):.2f} | {np.mean([r[c]['r_pool'][0] for r in G[t]]):.2f} |")
    return "\n".join(out)
def nc_cost():
    C = cols("0.004") + cols("0.012")
    out = ["| | " + " | ".join(f"{lab}, p={t.split('p=')[1]}" for t, c, lab in C) + " |", "|---|" + "---|" * len(C)]
    for key, lab in [("PEC full", "plain PEC"), ("PEC 1st order + ED", "1st-order PEC + detection"), ("ED", "error detection only"), ("ext two-point", "two-point extrapolation"),
                     ("u-amp exp fit", "ZNE, undetected amplified"), ("ZNE raw exp fit", "ZNE, no checks"), ("ZNE on ED exp fit", "ZNE on ED, all amplified")]:
        out.append(f"| {lab} | " + " | ".join(f"{np.median([r[c]['cost'][key] for r in G[t]]):.1f}" for t, c, _ in C) + " |")
    return "\n".join(out)

# ---------------- code appendix
FILES = ["edzne.py", "sweeps.py", "checks.py", "fig_subsets.py", "fig_scaling.py",
         "nonclifford/ising_ed.py", "nonclifford/study.py", "nonclifford/run_all.py", "nonclifford/report.py", "nonclifford/summary.py",
         "nonclifford/test_sim.py", "nonclifford/fig_nc.py", "nonclifford/fig_summary.py", "build_report.py"]
def code_appendix():
    out = []
    for f in FILES:
        src = open(P(f)).read().rstrip("\n")
        out.append(f"### `{f}`\n\n<details>\n<summary>{len(src.splitlines())} lines — click to expand</summary>\n\n```python\n{src}\n```\n\n</details>\n")
    return "\n".join(out)

if __name__ == "__main__":
    s = open(P("report_template.md")).read()
    for tag, val in [("{{CLIFFORD_TABLE}}", clifford_table()), ("{{NC_TABLE_LO}}", nc_table("0.004")), ("{{NC_TABLE_HI}}", nc_table("0.012")),
                     ("{{NC_RATIOS}}", nc_ratios()), ("{{NC_COST}}", nc_cost()), ("{{CODE_APPENDIX}}", code_appendix())]:
        assert tag in s, tag
        s = s.replace(tag, val)
    assert "{{" not in s.split("## Appendix: code")[0]
    open(P("REPORT.md"), "w").write(s)
    print(f"REPORT.md written: {len(s.splitlines())} lines")
