"""Figure for Stage 0 (fig_stage0.png) from stage0_results.json, plus the first/second-order split of the
last-round-only undetected rate computed from the detector error model."""
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from iceberg import *

R0 = json.load(open('stage0_results.json'))
C = dict(blue='#2a78d6', orange='#eb6834', aqua='#1baf7a', yellow='#eda100', gray='#52514e')
fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))

# (a) post-selected logical error rate
a = ax[0]; pp = np.logspace(-3.3, -1.4, 50)
for row, col, ls in ((r, c_, l) for r in R0['ii'] for c_, l in [((C['blue'] if r['prep'] == 'ft' else C['orange']), ('-' if r['R'] == 1 else '--'))] if r['R'] in (1, 8)):
    lab = f"{'FT' if row['prep']=='ft' else 'non-FT'} prep, R={row['R']}"
    a.plot(pp, row['A_Z'] * pp + row['B_Z'] * pp ** 2, ls, color=col, lw=2, label=lab + f" (A={row['A_Z']:.2f}, B={row['B_Z']:.0f})")
    ps = sorted(float(k) for k in row['sampled']); a.plot(ps, [row['sampled'][str(k) if str(k) in row['sampled'] else k]['pL'] for k in ps], 'o', color=col, ms=6)
a.plot(pp, pp, ':', color='#999', lw=1); a.text(1.2e-3, 1.6e-3, 'p_L = p', color='#777', fontsize=8)
a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('physical error p'); a.set_ylabel('post-selected logical error p_L (Z memory)')
a.set_title('(a) p_L = A p + B p²: lines exact (DEM), dots sampled', fontsize=9); a.legend(fontsize=7, frameon=False)

# (b) last-round-only amplification
b = ax[1]; Rs = [r['R'] for r in R0['iii']]
b.plot(Rs, [r['frac_part'] for r in R0['iii']], 'o-', color=C['blue'], lw=2, label='last round only (sampled)')
b.plot(Rs, [r['frac_full'] for r in R0['iii']], 's-', color=C['gray'], lw=2, label='full post-selection (sampled)')
first, second = [], []
for R in Rs:
    c, meta = memory_circuit("Z", R, Noise(3e-3), False, "ft", absolute=True)
    P, H, L = dem_mechanisms(c); lab = meta['labels']
    use = np.array([(k == "prep") or (k in "ZX" and r < 1) or k.startswith("abs") for k, r in lab]); Hp = H[:, use]
    first.append(P[(L.any(1)) & ~Hp.any(1)].sum())
    keys = {}
    for i in range(len(P)): keys.setdefault(Hp[i].tobytes(), []).append(i)
    s2 = 0.0
    for idx in keys.values():
        for x in range(len(idx)):
            for y in range(x + 1, len(idx)):
                if (L[idx[x]] ^ L[idx[y]]).any(): s2 += P[idx[x]] * P[idx[y]]
    second.append(s2)
b.plot(Rs, np.where(np.array(first) > 0, first, np.nan), '--', color=C['orange'], lw=2, label='1st order: unflagged hooks (∝R)')
b.plot(Rs, second, '--', color=C['aqua'], lw=2, label='2nd order: cancelling pairs (∝R²)')
b.set_xscale('log'); b.set_yscale('log'); b.set_xlabel('syndrome rounds R'); b.set_ylabel('undetected logical fraction')
b.set_title('(b) last-round-only amplification, p = 0.3%', fontsize=9); b.legend(fontsize=7, frameon=False)

# (c) geometric law in the hidden fault count k
cc = ax[2]
for g, col, mk in ((g, C['blue'] if g['stat'] == 'poisson' else C['orange'], 'o' if g['p'] == 0.01 else '^') for g in R0['iv']):
    k = np.arange(7); n = np.array(g['n_k'][:7]); E = np.array(g['E_k'][:7]); ok = n > 2000
    cc.errorbar(k[ok], E[ok], yerr=1 / np.sqrt(n[ok]), fmt=mk, color=col, ms=6, capsize=2, label=f"{g['stat']}, p={g['p']}")
g = R0['iv'][2]; kk = np.linspace(0, 6, 50); cc.plot(kk, g['O_u'] * g['rho'] ** kk, '-', color=C['gray'], lw=1.5, label=f"O_u ρ̄^k, ρ̄={g['rho']:.3f} (DEM)")
cc.set_yscale('log'); cc.set_xlabel('number of detected fault events k (hidden)'); cc.set_ylabel('E[o | k]')
cc.set_title('(c) geometric law, R = 8 rounds', fontsize=9); cc.legend(fontsize=7, frameon=False)

# (d) what an experiment sees: number of fired detectors
d = ax[3]
for g, col, mk in ((g, C['blue'] if g['p'] == 0.01 else C['yellow'], 'o') for g in R0['iv'] if g['stat'] == 'bernoulli'):
    f = np.arange(9); n = np.array(g['n_f'][:9]); E = np.array(g['E_f'][:9]); ok = n > 2000
    d.errorbar(f[ok], E[ok], yerr=1 / np.sqrt(n[ok]), fmt='o-', color=col, ms=5, capsize=2, label=f"p={g['p']}")
d.set_yscale('log'); d.set_xlabel('number of fired detectors |s| (observed)'); d.set_ylabel('E[o | |s|]')
d.set_title('(d) same data vs fired-detector count: not geometric', fontsize=9); d.legend(fontsize=7, frameon=False)
for x in ax: x.grid(alpha=0.25, lw=0.6); x.spines[['top', 'right']].set_visible(False)
fig.tight_layout(); fig.savefig('fig_stage0.png', dpi=140)
json.dump(dict(R=Rs, first=first, second=second), open('stage0_lastround_split.json', 'w'))
print('ok')
