"""Test of the cost-exponent model: ZNE two-point exponential pays (4 e^{2 Lam_O} + e^{4 Lam_O}) with Lam_O = ln(O_ideal/O_raw)
the observable's own decay; the twirled-check estimator pays 1/q_0^2 = e^{2 Lam_tot} with Lam_tot the full visible fault weight.
Uses the 1%, 0.5%, 0.2% runs; sigma^2 = 0.2 per shot for the 5-site mean of +-1 outcomes."""
import json, numpy as np
base = json.load(open("demo_results_2_4_6_8_10_12_14_16.json")); rep = json.load(open("repaired_sweep.json"))["depths"]; N, sig2 = 10000, 0.2
pts = []
for d, r in base["depths"].items():
    I = np.mean(r["ideal"][:5]); pts.append(dict(p=0.01, d=int(d), I=I, raw=r["stage1"]["raw"]["mx"], zne=r["stage1"]["ZNE exp (1,2)"]["mx_std"], cum2=r["stage1"]["ZNE cum2 (1,2,3)"]["mx_std"],
        rep=rep[d]["twirled stratified"]["mx_std"], q0r=rep[d]["species"][0], con=r["stage3"]["twirled stratified"]["mx_std"], q0c=r["stage3"]["species"][0]))
for p in (0.005, 0.002):
    for d, r in json.load(open(f"threshold_p{p}.json"))["depths"].items():
        pts.append(dict(p=p, d=int(d), I=r["ideal"], raw=r["raw"]["mx"], zne=r["ZNE exp (1,2)"]["mx_std"], cum2=r["ZNE cum2 (1,2,3)"]["mx_std"], rep=r["repaired full twirl"]["mx_std"], q0r=r["repaired full twirl"]["q0"], con=r["twirled stratified"]["mx_std"], q0c=r["twirled stratified"]["q0"]))
rows = ["| p | d | $\\Lambda_O$ | $\\Lambda_{tot}$ (rep.) | $\\Lambda_O/\\Lambda_{tot}$ | ZNE std pred / obs | repaired std pred / obs | constrained std pred / obs | ratio rep/ZNE pred / obs |", "|---|---|---|---|---|---|---|---|---|"]
for t in pts:
    LO = np.log(t["I"]/t["raw"]); Lt = -np.log(t["q0r"]); Lc = -np.log(t["q0c"])
    vz = (2*sig2/N)*(4*np.exp(2*LO)+np.exp(4*LO)); vr = (sig2/N)/t["q0r"]**2; vc = (sig2/N)/t["q0c"]**2
    rows.append(f"| {t['p']*100:g}% | {t['d']} | {LO:.2f} | {Lt:.2f} | {LO/Lt:.2f} | {np.sqrt(vz):.4f} / {t['zne']:.4f} | {np.sqrt(vr):.4f} / {t['rep']:.4f} | {np.sqrt(vc):.4f} / {t['con']:.4f} | {np.sqrt(vr/vz):.2f} / {t['rep']/t['zne']:.2f} |")
open("breakeven_check.md", "w").write("\n".join(rows) + "\n"); print("\n".join(rows))
