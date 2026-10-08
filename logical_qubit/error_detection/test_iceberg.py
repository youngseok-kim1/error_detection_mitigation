"""Stage 0 checks of iceberg.py.
1. every detector and observable is deterministic (stim raises otherwise) and a noiseless run never fires;
2. 1-FT: with perfect or FT preparation, no single fault flips a logical observable without firing a detector,
   for 1-4 rounds, with and without ancilla reset, in both memory bases;
3. non-FT preparation: list the first-order undetected logical faults (the A p term of the preparation)."""
import itertools, numpy as np
from iceberg import *

nz = Noise(1e-3)
print("1. determinism and noiseless runs")
for basis, R, reset, prep in itertools.product("ZX", (0, 1, 2, 3, 4), (False, True), ("perfect", "nft", "ft")):
    c, _ = memory_circuit(basis, R, nz, reset, prep)
    c.detector_error_model()                                         # raises on a non-deterministic detector/observable
    c0, _ = memory_circuit(basis, R, None, reset, prep)
    det, obs = c0.compile_detector_sampler().sample(200, separate_observables=True)
    assert not det.any() and not obs.any(), (basis, R, reset, prep)
print("   all variants deterministic; noiseless runs fire nothing")

print("2. single-fault check (perfect / FT preparation)")
for basis, R, reset, prep in itertools.product("ZX", (1, 2, 3, 4), (False, True), ("perfect", "ft")):
    bad = single_fault_report(memory_circuit(basis, R, nz, reset, prep)[0])
    assert not bad, (basis, R, reset, prep, bad)
print("   no first-order undetected logical error in 2 bases x 4 depths x reset/no-reset x 2 preparations")

print("3. non-FT preparation: first-order undetected logical mechanisms")
for basis in "ZX":
    for reset in (False, True):
        bad = single_fault_report(memory_circuit(basis, 2, nz, reset, "nft")[0])
        tot = sum(p for p, _ in bad)
        print(f"   basis {basis} reset={reset}: {len(bad)} mechanisms, total probability {tot:.2e} = {tot/nz.p:.2f} p  (observables {sorted(set(tuple(o) for _, o in bad))})")
print("all tests passed")
