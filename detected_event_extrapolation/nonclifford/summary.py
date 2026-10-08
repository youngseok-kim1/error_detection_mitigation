"""Compact summary table (rms |bias| of single-site <X_i>) from the .jsonl files"""
import json, numpy as np
from collections import defaultdict
recs = [json.loads(l) for fn in ("flagship.jsonl", "generic_lo.jsonl", "generic_hi.jsonl") for l in open(fn)]
G = defaultdict(list)
for r in recs: G[r["tag"]].append(r)
rows = [("raw", "unmitigated"), ("ED", "error detection only"), ("ext, coverage ratio", "extrapolated, coverage ratio"),
        ("ext, oracle 1st-order ratio", "extrapolated, exact 1st-order ratio"), ("ext, oracle ratio + pair term", "  ... + in-situ pair term"),
        ("ext, site-pooled training ratio", "extrapolated, Clifford-trained ratio"), ("CDR on ED (pooled slope)", "Clifford regression on ED data"),
        ("PEC 1st order (exact model)", "1st-order PEC, exact model"), ("PEC 1st order (model 30% high)", "1st-order PEC, model 30% high"),
        ("u-amp ZNE, exp fit (G=1,2,3)", "ZNE, undetected faults amplified"), ("ZNE on ED, exp fit (G=1,2,3)", "ZNE on ED, all noise amplified"),
        ("ZNE raw, exp fit (G=1,2,3)", "ZNE, no checks")]
cols = [("flag d=6 p=0.004", "anc"), ("flag d=6 p=0.004", "parity"), ("generic p=0.004", "anc"), ("generic p=0.004", "parity"),
        ("flag d=6 p=0.012", "anc"), ("flag d=6 p=0.012", "parity"), ("generic p=0.012", "anc"), ("generic p=0.012", "parity")]
def rms(tag, cfg, key):
    n = G[tag][0]["n"]; B = np.array([np.array(r[cfg]["est"][key])[:n] - np.array(r["ideal"])[:n] for r in G[tag]]); return np.sqrt((B ** 2).mean())
print("rms |bias| of <X_i> over sites (and instances)\n")
print(f"{'':40s}" + "".join(f"{t.replace('flag ','paper ').replace(' p=',' p='):>20s}" for t, _ in cols))
print(f"{'':40s}" + "".join(f"{('anc only' if c=='anc' else 'anc+parity'):>20s}" for _, c in cols))
for key, lab in [("Lam", "total fault rate Lambda"), ("coverage", "coverage"), ("alpha", "acceptance")]:
    print(f"{lab:40s}" + "".join(f"{np.mean([(r[c]['Ld']+r[c]['Lu']) if key=='Lam' else r[c][key] for r in G[t]]):20.3f}" for t, c in cols))
print(f"{'mean |ideal <X_i>|':40s}" + "".join(f"{np.mean(np.abs([r['ideal'][:r['n']] for r in G[t]])):20.3f}" for t, c in cols))
for key, lab in rows: print(f"{lab:40s}" + "".join(f"{rms(t, c, key):20.4f}" for t, c in cols))
print("\nsampling cost, N x Var for one +-1 observable (median over instances; training shots not included)")
for key in ("PEC full", "PEC 1st order + ED", "ED", "ext two-point", "u-amp exp fit", "ZNE raw exp fit", "ZNE on ED exp fit"):
    print(f"{key:40s}" + "".join(f"{np.median([r[c]['cost'][key] for r in G[t]]):20.1f}" for t, c in cols))
print("\nfitted ratios: coverage ratio | exact 1st-order (mean over sites) | Clifford-trained (pooled)")
for t, c in cols: print(f"  {t:22s} {c:7s} {np.mean([r[c]['r_cov'] for r in G[t]]):.3f} | {np.mean([np.mean(r[c]['r1'][:r['n']]) for r in G[t]]):.3f} | {np.mean([r[c]['r_pool'][0] for r in G[t]]):.3f}")
a = [abs(v - r[c]["alpha"]) for r in recs for c in ("anc", "parity") for v in r[c]["alpha_uamp"].values()]
print(f"\nacceptance change when only undetected faults are amplified: max {max(a):.1e}")
