"""Exact OPT for one-round 2-set agreement, three processes, three values, independent links (second route, no solver).

On an input with three distinct values a view that hears one message already determines the input, so the free views of
different inputs are disjoint and every input can be treated alone; by relabelling all six are the same problem as (0,1,2).
On it a rule picks, for each process, a source: hear nothing -> self; hear only j -> self or j; hear only k -> self or k;
hear both -> self, j or k. That is 12 choices per process and 1728 rules. A pattern fails when the three sources are
distinct. For each rule we count failing patterns by number of delivered links (exact integers), drop every count vector
that is coefficientwise at least another (it is never strictly better for any p), and compute the lower envelope of the
rest on [0,1] with exact real roots (sympy).

usage: python env650.py
"""
import itertools
import sympy as sp

P = (0, 1, 2)
LINKS = [(s, r) for r in P for s in P if s != r]


def choices(i):
    j, k = [x for x in P if x != i]
    for a in (i, j):              # hear only j
        for b in (i, k):          # hear only k
            for c in (i, j, k):   # hear both
                yield {(): i, (j,): a, (k,): b, (j, k): c}


def counts(rule):
    cnt = [0] * 7
    for mask in range(64):
        got = {LINKS[t] for t in range(6) if mask >> t & 1}
        src = []
        for i in P:
            heard = tuple(s for s in P if s != i and (s, i) in got)
            src.append(rule[i][heard])
        if len(set(src)) == 3:
            cnt[bin(mask).count("1")] += 1
    return tuple(cnt)


def main():
    per = [list(choices(i)) for i in P]
    vecs = {}
    for r in itertools.product(*per):
        vecs.setdefault(counts(r), r)
    print(f"rules 1728, distinct count vectors {len(vecs)}")
    keep = [v for v in vecs if not any(w != v and all(a <= b for a, b in zip(w, v)) for w in vecs)]
    print(f"not coefficientwise dominated: {len(keep)}")
    p = sp.symbols("p")
    poly = {v: sp.expand(sum(c * p ** k * (1 - p) ** (6 - k) for k, c in enumerate(v))) for v in keep}
    # breakpoints: all real roots in (0,1) of pairwise differences, then pick the minimiser on each piece
    pts = {sp.Integer(0), sp.Integer(1)}
    for a, b in itertools.combinations(keep, 2):
        d = sp.Poly(poly[a] - poly[b], p)
        if d.is_zero:
            continue
        for rt in d.real_roots():
            if 0 < rt < 1:
                pts.add(rt)
    pts = sorted(pts, key=lambda x: float(x))
    pieces = []
    for lo, hi in zip(pts, pts[1:]):
        mid = (lo + hi) / 2
        best = min(keep, key=lambda v: float(poly[v].subs(p, mid).evalf(50)))
        if pieces and pieces[-1][2] == best:
            pieces[-1][1] = hi
        else:
            pieces.append([lo, hi, best])
    print("lower envelope (OPT on [0,1]):")
    for lo, hi, v in pieces:
        print(f"  [{float(lo):.6f}, {float(hi):.6f}]  counts {v}  error {sp.factor(poly[v])}  ends {lo} .. {hi}")
        r = vecs[v]
        for i in P:
            print(f"      process {i}: " + ", ".join(f"hears {list(h) if h else 'none'} -> {r[i][h]}" for h in r[i]))
    grid = ["1/100", "1/20", "1/10", "1/5", "1/4", "1/3", "2/5", "1/2", "3/5", "2/3", "3/4", "4/5", "9/10", "19/20", "99/100"]
    print("envelope at the CP-SAT grid points:")
    for g in grid:
        x = sp.Rational(g)
        val = min(poly[v].subs(p, x) for v in keep)
        print(f"  p={g} OPT={val}")


if __name__ == "__main__":
    main()
