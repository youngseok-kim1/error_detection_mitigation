"""Tabulate biases from the .jsonl outputs of run_all.py"""
import sys, json, numpy as np
from collections import defaultdict
def load(fn): return [json.loads(l) for l in open(fn)]
def bias(rec, cfg, name):
    return np.array(rec[cfg]["est"][name]) - np.array(rec["ideal"])
def show(recs, title):
    n = recs[0]["n"]; names = list(recs[0]["anc"]["est"])
    for cfg, lab in (("anc", "ancilla checks only"), ("parity", "ancilla checks + parity")):
        c = [r[cfg] for r in recs]
        print(f"\n### {title} | {lab} | {len(recs)} instance(s)")
        print(f"    Lam={np.mean([x['Ld']+x['Lu'] for x in c]):.3f}  coverage={np.mean([x['coverage'] for x in c]):.3f}  alpha={np.mean([x['alpha'] for x in c]):.3f}"
              f"  r_cov={c[0]['r_cov']:.3f}  r1(X_i) mean={np.mean([x['r1'][:n] for x in c]):.3f} [{np.min([x['r1'][:n] for x in c]):.2f},{np.max([x['r1'][:n] for x in c]):.2f}]"
              f"  r_train(pooled)={np.mean([x['r_pool'][0] for x in c]):.3f}")
        print(f"    |ideal <X_i>| mean={np.mean(np.abs([r['ideal'][:n] for r in recs])):.3f}")
        print("    cost (N x Var, central site): " + "  ".join(f"{k}={np.median([x['cost'][k] for x in c]):.1f}" for k in c[0]["cost"]))
        print(f"    {'estimator':38s} {'X_i rms':>9s} {'X_i max':>9s} {'m_x mean':>10s} {'m_x rms':>9s} {'XX rms':>9s}")
        for nm in names:
            B = np.array([bias(r, cfg, nm) for r in recs]); bx, bm, bxx = B[:, :n], B[:, n], B[:, n + 1:2 * n + 1]
            print(f"    {nm:38s} {np.sqrt((bx**2).mean()):9.5f} {np.abs(bx).max():9.5f} {bm.mean():+10.5f} {np.sqrt((bm**2).mean()):9.5f} {np.sqrt((bxx**2).mean()):9.5f}")
if __name__ == "__main__":
    for fn in sys.argv[1:]:
        recs = load(fn); groups = defaultdict(list)
        for r in recs: groups[r["tag"]].append(r)
        for tag, g in groups.items(): show(g, tag)
