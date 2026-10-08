"""Flag statistics for the heavy-hex design: bond (i,i+1) mediated by m (gadget CX(i->m) CX(i+1->m) Rz(m) CX(i+1->m) CX(i->m)),
mediator measured (fixed check, image type t on data and on m), peripheral checks p_i, p_{i+1} (single-wire, image type t_i, t_{i+1}, enabled w.p. b).
Compilation twirl: image type t drawn uniformly from {X,Y,Z}, the same for the mediator's data components and the peripheral image on that qubit
(the compilation fixes the frame of the data qubit at the noise location). Noise: 2-qubit depolarizing on (q,m) after each CX."""
import itertools, numpy as np
AX = 'XYZ'
def anti(a, b): return a != 'I' and b != 'I' and a != b
P2 = [a + b for a in 'IXYZ' for b in 'IXYZ' if a + b != 'II']
# locations: which data qubit the CX touches and the mediator image's data support at that time
LOCS = {"after CX(i->m) #1": ('i', ['i']), "after CX(i+1->m) #2": ('j', ['i', 'j']), "after CX(i+1->m) #3": ('j', ['i']), "after CX(i->m) #4": ('i', [])}
def stats(b, use_med=True):
    out = {}
    for name, (dq, msup) in LOCS.items():
        dist = {}
        for pp in P2:
            fd, fm = pp[0], pp[1]                   # fault on data qubit dq and on mediator
            fl = np.zeros(4)                         # distribution of number of flags 0..3
            for ti, tj, tm in itertools.product(AX, AX, AX):
                t = {'i': ti, 'j': tj}
                # mediator image type on data follows the data frame; on m it is tm (mediator type fixed by compilation; take Z-type vs X-type -> uniform ok)
                med = (sum(anti(fd if q == dq else 'I', t[q]) for q in msup) + anti(fm, tm)) % 2 if use_med else 0
                per = anti(fd, t[dq])                # the peripheral check on the faulty data qubit
                for e in (0, 1):
                    pe = b if e else 1 - b
                    fl[med + (e and per)] += pe / 27
            dist[pp] = fl
        out[name] = dist
    return out
for b, um, label in ((0, True, "mediator checks only (anchor-style, type-twirled)"), (0.75, True, "mediator + peripheral b=3/4"), (0.75, False, "peripheral only b=3/4 (mediator flags ignored)")):
    print(f"\n== {label} ==")
    for name, dist in stats(b, um).items():
        p1 = {pp: round(1 - fl[0], 3) for pp, fl in dist.items()}
        vals = sorted(set(p1.values()))
        print(f"  {name}: P(>=1 flag) values {vals}; mean {np.mean(list(p1.values())):.3f}; invisible: {[pp for pp, v in p1.items() if v == 0]}")
