"""Optimal deterministic one-round consensus for n processes, q input values,
independent links, as a mixed integer program at a rational p.

Variables: d[i,view,v] in {0,1}, one value per view, v ranging over the values the
view contains (strong validity: a deterministic rule must decide a value it has
seen). y[x,pat,v] in [0,1] with y <= d[i, view_i(x,pat), v] for every i: y is the
indicator that all processes decide v. Minimise t subject to, for every
non-constant input x, t >= 1 - sum_pat w(pat) sum_v y[x,pat,v].

The optimal rule found is then re-evaluated exactly in rationals (Fractions), and
its per-input disagreement polynomials are printed as counts by number of
delivered links. The MILP value is a floating claim; the printed exact value is
the value of an explicit algorithm, hence an exact upper bound.

usage: python milp3.py n q a/b [a/b ...]
"""
import sys
import itertools
from fractions import Fraction
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix


def setup(n, q):
    others = {i: [j for j in range(n) if j != i] for i in range(n)}
    links = [(j, i) for i in range(n) for j in others[i]]
    lidx = {l: k for k, l in enumerate(links)}
    inputs = [x for x in itertools.product(range(q), repeat=n) if len(set(x)) > 1]
    return others, links, lidx, inputs


def view(i, x, pat, others, lidx):
    return (x[i],) + tuple(x[j] if (pat >> lidx[(j, i)]) & 1 else None for j in others[i])


def seen(v):
    return sorted(set(s for s in v if s is not None))


def solve(n, q, p, time_limit=3600):
    others, links, lidx, inputs = setup(n, q)
    L = len(links)
    pats = list(range(1 << L))
    w = [p ** bin(pt).count("1") * (1 - p) ** (L - bin(pt).count("1")) for pt in pats]
    # decision variables only for views with >= 2 seen values
    dvar = {}
    nv = 0
    for x in inputs:
        for pt in pats:
            for i in range(n):
                vw = view(i, x, pt, others, lidx)
                sv = seen(vw)
                if len(sv) >= 2:
                    for v in sv:
                        if (i, vw, v) not in dvar:
                            dvar[(i, vw, v)] = nv
                            nv += 1
    nd = nv
    yvar = {}
    for xi, x in enumerate(inputs):
        for pt in pats:
            vws = [view(i, x, pt, others, lidx) for i in range(n)]
            common = set(range(q))
            for vw in vws:
                common &= set(seen(vw))
            for v in sorted(common):
                yvar[(xi, pt, v)] = nv
                nv += 1
    t = nv
    nv += 1
    rows = []  # (coeff dict, lo, hi)
    # one value per free view
    views = {}
    for (i, vw, v), k in dvar.items():
        views.setdefault((i, vw), []).append(k)
    for ks in views.values():
        rows.append(({k: 1.0 for k in ks}, 1.0, 1.0))
    # y <= d (d fixed to 1 for forced views, so no constraint there)
    for (xi, pt, v), k in yvar.items():
        x = inputs[xi]
        for i in range(n):
            vw = view(i, x, pt, others, lidx)
            if (i, vw, v) in dvar:
                rows.append(({k: 1.0, dvar[(i, vw, v)]: -1.0}, -np.inf, 0.0))
    # t + sum w y >= 1
    for xi in range(len(inputs)):
        c = {t: 1.0}
        for pt in pats:
            for v in range(q):
                k = yvar.get((xi, pt, v))
                if k is not None:
                    c[k] = c.get(k, 0.0) + float(w[pt])
        rows.append((c, 1.0, np.inf))
    A = lil_matrix((len(rows), nv))
    lo = np.empty(len(rows))
    hi = np.empty(len(rows))
    for r, (c, a, b) in enumerate(rows):
        for k, val in c.items():
            A[r, k] = val
        lo[r], hi[r] = a, b
    cvec = np.zeros(nv)
    cvec[t] = 1.0
    integrality = np.zeros(nv)
    integrality[:nd] = 1
    res = milp(cvec, constraints=LinearConstraint(A.tocsr(), lo, hi), integrality=integrality,
               bounds=Bounds(np.zeros(nv), np.ones(nv)),
               options={"time_limit": time_limit, "mip_rel_gap": 0.0, "mip_feasibility_tolerance": 1e-9})
    rule = {}
    for (i, vw), ks in views.items():
        best = max(ks, key=lambda k: res.x[k])
        for (ii, vv, v), k in dvar.items():
            if k == best:
                rule[(i, vw)] = v
    return res, rule, inputs, (others, links, lidx)


def exact_eval(n, q, rule, inputs, geo):
    others, links, lidx = geo
    L = len(links)
    out = []
    for x in inputs:
        cnt = [0] * (L + 1)
        for pt in range(1 << L):
            ds = []
            for i in range(n):
                vw = view(i, x, pt, others, lidx)
                sv = seen(vw)
                ds.append(sv[0] if len(sv) == 1 else rule[(i, vw)])
            if len(set(ds)) > 1:
                cnt[bin(pt).count("1")] += 1
        out.append((x, cnt))
    return out


def poly_at(cnt, p):
    L = len(cnt) - 1
    return sum(c * p ** k * (1 - p) ** (L - k) for k, c in enumerate(cnt))


def main():
    n, q = int(sys.argv[1]), int(sys.argv[2])
    for s in sys.argv[3:]:
        p = Fraction(s)
        res, rule, inputs, geo = solve(n, q, float(p))
        ev = exact_eval(n, q, rule, inputs, geo)
        vals = [(poly_at(c, p), x, c) for x, c in ev]
        worst = max(v for v, _, _ in vals)
        print(f"n={n} q={q} p={p} milp_status={res.status} milp_obj={res.fun:.9f} "
              f"mip_dual_bound={getattr(res, 'mip_dual_bound', float('nan')):.9f} "
              f"exact_rule_value={worst} ({float(worst):.9f})")
        seen_c = {}
        for v, x, c in vals:
            seen_c.setdefault(tuple(c), []).append("".join(map(str, x)))
        for c, xs in sorted(seen_c.items(), key=lambda kv: -float(poly_at(list(kv[0]), p))):
            print(f"   counts {list(c)} value {float(poly_at(list(c), p)):.6f} inputs {' '.join(xs)}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
