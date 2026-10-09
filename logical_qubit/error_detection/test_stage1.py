"""Stage 1 checks.
1. noiseless encoded and unencoded trajectories reproduce the exact logical Trotter values;
2. noiseless encoded runs fire no detector and no terminal parity, for every syndrome schedule;
3. skeleton claim: with noise, the detector statistics of the non-Clifford trajectories equal those of the stim
   skeleton (rotations dropped), detector by detector and for pairs of detectors;
4. fault-count estimator: k_hat = 0 iff no detector fired; every single DEM mechanism has k_hat = 1."""
import numpy as np
from circuits import *
from sim_mc import sample
from estimators import *

print("1/2. noiseless runs")
for steps, every in ((1, 1), (3, 1), (4, 2), (5, 'end'), (6, None)):
    c = encoded_mfim(steps, every); rec = sample(c, Params(0, 0, 0), 100_000, seed=steps)
    V, par, o, _ = unpack(rec, c); I = ideal_logical(steps)
    for name in ('Z0', 'Z1', 'ZZ'):
        val = unpack(rec, c, name)[2].mean(); se = np.sqrt(max(1 - I[name] ** 2, 1e-6) / 1e5)
        assert abs(val - I[name]) < 5 * se + 1e-9, (steps, every, name, val, I[name])
    assert not V.any() and not par.any(), (steps, every)
    u = unencoded_mfim(steps); ru = sample(u, Params(0, 0, 0), 100_000, seed=steps)
    uv = unencoded_obs(ru, 'Z0').mean(); assert abs(uv - I['Z0']) < 5 * np.sqrt((1 - I['Z0'] ** 2 + 1e-6) / 1e5)
    print(f"   steps={steps} every={every}: Zbar0 {val if name=='ZZ' else 0:.0f}.. ok (ideal Z0 {I['Z0']:+.4f}, ZZ {I['ZZ']:+.4f}); detectors silent; unencoded ok")

print("3. skeleton vs trajectories (p = 1%)")
prm = Params(1e-2)
for steps, every in ((3, 1), (4, 2)):
    c = encoded_mfim(steps, every); N = 400_000
    V = unpack(sample(c, prm, N, seed=11), c)[0].astype(float)
    sk = skeleton(c, prm)[0]; Vs = sk.compile_detector_sampler().sample(N).astype(float)
    m1, m2 = V.mean(0), Vs.mean(0); se = np.sqrt(m2 * (1 - m2) * 2 / N) + 1e-9
    C1, C2 = (V.T @ V) / N, (Vs.T @ Vs) / N; seC = np.sqrt(C2 * (1 - C2) * 2 / N) + 1e-9
    z1 = np.abs(m1 - m2) / se; z2 = np.abs(C1 - C2) / seC
    print(f"   steps={steps} every={every}: {V.shape[1]} detectors, max |z| single {z1.max():.1f}, pairs {z2.max():.1f}")
    assert z1.max() < 5 and z2.max() < 5.5

print("4. k_hat")
c = encoded_mfim(4, 1); kh = KHat(c, Params(3e-3))
dem = skeleton(c, Params(3e-3))[0].detector_error_model(); K = dem.num_detectors
for inst in dem.flattened():
    if inst.type != 'error': continue
    v = np.zeros((1, K), np.uint8)
    for t in inst.targets_copy():
        if t.is_relative_detector_id(): v[0, t.val] ^= 1
    assert kh(v)[0] == (1 if v.any() else 0)
print(f"   {kh.n_mech} distinct detector signatures, {K} detectors; every single mechanism has k_hat = 1")
print("all tests passed")

print("5. exact density matrix vs trajectories (p = 1%)")
from dm_exact import values
from circuits import classify_faults
for steps, every in ((2, 1), (3, 'end')):
    c = encoded_mfim(steps, every); prm = Params(1e-2); N = 600_000
    und, harm = classify_faults(c)
    for zero, lab in ((None, 'noisy'), (harm, 'PEC (harmful undetected removed)')):
        ex = values(c, prm, ('Z0', 'ZZ'), 'ps', zero); ex_all = values(c, prm, ('Z0',), 'all', zero)
        rec = sample(c, prm, N, seed=5, zero=zero)
        for name in ('Z0', 'ZZ'):
            V, par, o, _ = unpack(rec, c, name)
            acc = (par == 0) & ~V.any(1); mc = o[acc].mean(); se = o[acc].std() / np.sqrt(acc.sum())
            assert abs(mc - ex[name]['o_par']) < 5 * se, (steps, every, lab, name, mc, ex[name]['o_par'])
            assert abs(acc.mean() - ex[name]['alpha_par']) < 5 * np.sqrt(acc.mean() * (1 - acc.mean()) / N)
        V, par, o, _ = unpack(rec, c, 'Z0'); acc = (par == 0) & ~V.any(1); mc = o[acc].mean()
        mc_raw = o.mean()
        assert abs(mc_raw - ex_all['Z0']['o_all']) < 5 / np.sqrt(N)
        print(f"   steps={steps} every={every} {lab}: PS <Zbar0> exact {ex['Z0']['o_par']:+.4f} vs MC {mc:+.4f}; acceptance {ex['Z0']['alpha_par']:.4f} vs {acc.mean():.4f}")
    u = unencoded_mfim(steps); eu = values(u, prm, ('Z0', 'ZZ'), 'all'); ru = sample(u, prm, N, seed=6)
    assert abs(unencoded_obs(ru, 'Z0').mean() - eu['Z0']['o_all']) < 5 / np.sqrt(N)
print("   exact DM agrees with trajectories (PS, acceptance, raw, PEC, unencoded)")
print("all tests passed (incl. 5)")
