"""Stage 0 experiments on one [[4,2,2]] block (memory, Clifford, stim). Output: stage0_results.json, fig_stage0.png.

(ii)  post-selected logical error rate p_L(p) for FT vs non-FT preparation and R syndrome rounds: exact first-order term A
      from the single-fault enumeration, second-order term B from all pairs of DEM mechanisms with cancelling syndrome
      and non-trivial logical effect, checked against direct stim sampling at larger p.
(iii) "last round only" amplification (Froland et al. App. B, Fig. 12c): R rounds, post-select only on the first round
      and on the last round + terminal readout; the undetected fraction should grow ~R^2 in an intermediate window,
      against ~R for full post-selection.
(iv)  geometric law E[o | k] = O_u rho^k (theory/qed_extrapolation.md, Theorem 1), with k the number of detected fault
      events, tested by sampling the DEM mechanisms directly with Poisson (Pauli-Lindblad) and Bernoulli (circuit-level)
      statistics; then what an experiment sees: E[o | number of fired detectors]."""
import json, numpy as np, stim, itertools, time
from iceberg import *

import sys, os
rng = np.random.default_rng(7)
PARTS = sys.argv[1].split(',') if len(sys.argv) > 1 else ['ii', 'iii', 'iv']
OUT = json.load(open('stage0_results.json')) if os.path.exists('stage0_results.json') else {}

# ---------------------------------------------------------------- (ii) p_L(p) = A p + B p^2
def AB(basis, R, prep, reset=False):
    p0 = 1e-3
    P, H, L = dem_mechanisms(memory_circuit(basis, R, Noise(p0), reset, prep)[0])
    lo = L.any(1); det = H.any(1)
    A = P[lo & ~det].sum() / p0
    # second order: pairs (a, b) with H_a = H_b (cancelling syndrome) and L_a xor L_b != 0
    keyH = [h.tobytes() for h in H]; B = 0.0
    groups = {}
    for i, k in enumerate(keyH): groups.setdefault(k, []).append(i)
    for k, idx in groups.items():
        if not det[idx[0]]: continue                                   # pairs of undetected first-order ones are O(p^2 * A), skip
        for a, b in itertools.combinations(idx, 2):
            if (L[a] ^ L[b]).any(): B += P[a] * P[b]
    return A, B / p0 ** 2

def sample_pL(basis, R, prep, p, shots, reset=False):
    c, _ = memory_circuit(basis, R, Noise(p), reset, prep)
    det, obs = c.compile_detector_sampler().sample(shots, separate_observables=True)
    acc = ~det.any(1); return float(obs[acc].any(1).mean()), float(acc.mean())

t0 = time.time(); res = []
for prep in (("ft", "nft") if 'ii' in PARTS else ()):
    for R in (1, 2, 4, 8):
        A, B = AB("Z", R, prep); AX_, BX = AB("X", R, prep)
        row = dict(prep=prep, R=R, A_Z=A, B_Z=B, A_X=AX_, B_X=BX, sampled={})
        for p in (3e-3, 1e-2, 3e-2):
            pl, al = sample_pL("Z", R, prep, p, 4_000_000)
            row["sampled"][p] = dict(pL=pl, alpha=al, pred=A * p + B * p * p)
        res.append(row)
        print(f"(ii) prep={prep} R={R}: Z-memory A={A:.3f} B={B:.1f} | X-memory A={AX_:.3f} B={BX:.1f} | sampled pL(Z) at p=1e-2: "
              f"{row['sampled'][1e-2]['pL']:.2e} vs Ap+Bp^2 {row['sampled'][1e-2]['pred']:.2e}, acceptance {row['sampled'][1e-2]['alpha']:.3f}", flush=True)
if "ii" in PARTS: OUT["ii"] = res

# ---------------------------------------------------------------- (iii) last-round-only amplification
def last_round_only(R, p, shots, keep_first=1, reset=False):
    c, meta = memory_circuit("Z", R, Noise(p), reset, "ft", absolute=True)
    det, obs = c.compile_detector_sampler().sample(shots, separate_observables=True)
    lab = meta["labels"]
    # [R1] App. B: a constant number of initial syndromes + the absolute value of the last round and of the terminal parity
    use_part = np.array([(k == "prep") or (k in "ZX" and r < keep_first) or k.startswith("abs") for k, r in lab])
    lab_std = np.array([not k.startswith("abs") for k, r in lab]); det_full = det[:, lab_std]
    acc_part = ~det[:, use_part].any(1); acc_full = ~det_full.any(1)
    return dict(R=R, frac_part=float(obs[acc_part].any(1).mean()), acc_part=float(acc_part.mean()),
                frac_full=float(obs[acc_full].any(1).mean()) if acc_full.any() else float('nan'), acc_full=float(acc_full.mean()))
res3 = []
for R in ((2, 4, 6, 8, 12, 16, 20, 30, 40) if 'iii' in PARTS else ()):
    res3.append(last_round_only(R, 3e-3, 2_000_000)); r = res3[-1]
    print(f"(iii) R={R}: last-round-only undetected fraction {r['frac_part']:.2e} (acc {r['acc_part']:.3f}) | full PS {r['frac_full']:.2e} (acc {r['acc_full']:.3f})", flush=True)
if "iii" in PARTS: OUT["iii"] = res3

# ---------------------------------------------------------------- (iv) geometric law
def geometric(R, p, shots, stat, reset=False, batch=200_000):
    c, meta = memory_circuit("Z", R, Noise(p), reset, "ft", obs=(0,))
    P, H, L = dem_mechanisms(c); L = L[:, 0]; det = H.any(1)
    lam = -0.5 * np.log(1 - 2 * P)
    O_u = float(np.exp(-2 * lam[~det & (L == 1)].sum()))
    rho = float((lam[det] * (1 - 2 * L[det])).sum() / lam[det].sum()); Lam_d = float(lam[det].sum())
    kmax = 12; sk = np.zeros(kmax + 1); nk = np.zeros(kmax + 1); sf = np.zeros(4 * R + 8); nf = np.zeros(4 * R + 8)
    s0 = n0 = 0.0
    for _ in range(shots // batch):
        n = rng.poisson(lam, size=(batch, len(lam))) if stat == "poisson" else (rng.random((batch, len(lam))) < P).astype(np.int64)
        k = (n[:, det]).sum(1); o = 1 - 2 * ((n @ L) % 2); s = (n @ H) % 2; f = s.sum(1)
        kk = np.minimum(k, kmax); np.add.at(sk, kk, o); np.add.at(nk, kk, 1)
        np.add.at(sf, f, o); np.add.at(nf, f, 1)
        z = f == 0; s0 += o[z].sum(); n0 += z.sum()
    Ek = sk / np.maximum(nk, 1); Ef = sf / np.maximum(nf, 1)
    return dict(R=R, p=p, stat=stat, O_u=O_u, rho=rho, Lam_d=Lam_d, alpha0_pred=float(np.exp(-Lam_d)),
                k=list(range(kmax + 1)), E_k=Ek.tolist(), n_k=nk.tolist(), pred_k=[O_u * rho ** k for k in range(kmax + 1)],
                E_f=Ef[:12].tolist(), n_f=nf[:12].tolist(), E_s0=s0 / n0, frac_s0=n0 / shots)
res4 = []
for p in ((3e-3, 1e-2) if 'iv' in PARTS else ()):
    for stat in ("poisson", "bernoulli"):
        g = geometric(8, p, 2_000_000, stat); res4.append(g)
        dev = [(g['E_k'][k] - g['pred_k'][k]) for k in range(6)]
        se = [1 / np.sqrt(max(g['n_k'][k], 1)) for k in range(6)]
        print(f"(iv) p={p} {stat}: Lambda_d={g['Lam_d']:.3f} rho={g['rho']:.4f} O_u={g['O_u']:.5f} | E[o|k]-pred, k=0..5: "
              + " ".join(f"{d:+.4f}({s:.4f})" for d, s in zip(dev, se)) + f" | E[o|s=0]={g['E_s0']:.5f}", flush=True)
if "iv" in PARTS: OUT["iv"] = res4
json.dump(OUT, open("stage0_results.json", "w"), indent=1, default=float)
print(f"done in {time.time()-t0:.0f}s")
