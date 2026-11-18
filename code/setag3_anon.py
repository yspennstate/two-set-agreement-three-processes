"""ANONYMOUS VARIANT: decisions tied across processes by (own value, multiset of heard values). Exact optimum of deterministic one-round k-set agreement (n processes, q input values, every directed link
delivers independently with probability p) at rational p = a/b, by CP-SAT in integer arithmetic.

A process decides a value it has seen (its own or one heard). An execution fails when more than k distinct values
are decided. The error of a rule is the largest failure probability over inputs; inputs with at most k distinct
values cannot fail, so only those with more than k distinct values enter. At p = a/b a pattern with j delivered
links out of L has weight a^j (b-a)^(L-j), total b^L. For each such input x:  t >= sum_pat w(pat) [fail].
Fail is modelled from below: ok_pat (at most k distinct decisions) needs, for n = 3 and k = 2, some pair i<j and
value v with both deciding v. CP-SAT minimises t; OPTIMAL is a proof; controls t <= OPT-1 INFEASIBLE, t <= OPT
FEASIBLE; the returned rule is re-evaluated by direct simulation (evaluate()), with no shared model code.

usage: python setag3.py n q k a/b [a/b ...]      (the pair encoding is written for k = n - 1)
"""
import itertools
import json
import os
import sys
from fractions import Fraction
from ortools.sat.python import cp_model
import milp3 as M


def setup(n, q, k):
    others, links, lidx, inputs = M.setup(n, q)
    inputs = [x for x in inputs if len(set(x)) > k]
    return others, links, lidx, inputs


def model(n, q, k, p, cap=None):
    assert k == n - 1, "pair encoding: fail = all n decisions distinct"
    others, links, lidx, inputs = setup(n, q, k)
    L = len(links)
    a, b = p.numerator, p.denominator
    pats = range(1 << L)
    w = [a ** bin(pt).count("1") * (b - a) ** (L - bin(pt).count("1")) for pt in pats]
    m = cp_model.CpModel()
    d, views, shared, tied = {}, {}, {}, set()

    def dec(i, vw, v):
        sv = M.seen(vw)
        if len(sv) == 1:
            return 1 if sv[0] == v else 0
        return d[(i, vw, v)]

    for x in inputs:
        for pt in pats:
            for i in range(n):
                vw = M.view(i, x, pt, others, lidx)
                sv = M.seen(vw)
                if len(sv) >= 2 and (i, vw) not in views:
                    # anonymous: the decision depends only on (own value, multiset of heard values), the same at every process
                    key = (vw[0], tuple(sorted(h for h in vw[1:] if h is not None)))
                    views[(i, vw)] = sv
                    for v in sv:
                        if (key, v) not in shared:
                            shared[(key, v)] = m.NewBoolVar("")
                        d[(i, vw, v)] = shared[(key, v)]
                    if key not in tied:
                        tied.add(key)
                        m.AddExactlyOne([shared[(key, v)] for v in sv])
    t = m.NewIntVar(0, b ** L, "t")
    for x in inputs:
        okterms = []
        for pt in pats:
            vws = [M.view(i, x, pt, others, lidx) for i in range(n)]
            z = []
            for i, j in itertools.combinations(range(n), 2):
                for v in sorted(set(M.seen(vws[i])) & set(M.seen(vws[j]))):
                    di, dj = dec(i, vws[i], v), dec(j, vws[j], v)
                    if isinstance(di, int) and di == 0 or isinstance(dj, int) and dj == 0:
                        continue
                    y = m.NewBoolVar("")
                    for dd in (di, dj):
                        if not isinstance(dd, int):
                            m.AddImplication(y, dd)
                    z.append(y)
            if not z:
                continue
            ok = m.NewBoolVar("")
            m.AddBoolOr(z).OnlyEnforceIf(ok)  # ok -> some pair agrees
            okterms.append((w[pt], ok))
        m.Add(t + sum(c * y for c, y in okterms) >= b ** L)
    if cap is None:
        m.Minimize(t)
    else:
        m.Add(t <= cap)
    return m, d, views, t, inputs, (others, links, lidx), b ** L


def run(n, q, k, p, cap=None, workers=2, tl=36000):
    m, d, views, t, inputs, geo, W = model(n, q, k, p, cap)
    s = cp_model.CpSolver()
    s.parameters.num_search_workers = workers
    s.parameters.max_time_in_seconds = tl
    st = s.Solve(m)
    rule = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        rule = {key: v for key, sv in views.items() for v in sv if s.Value(d[(key[0], key[1], v)])}
    return s.StatusName(st), (s.Value(t) if rule is not None else None), rule, inputs, geo, W


def evaluate(n, q, k, rule, p):
    """Direct simulation: for each input with > k values, failure counts by number of delivered links, and the
    worst failure probability at p."""
    procs = range(n)
    links = [(s, r) for r in procs for s in procs if s != r]
    worst, rows = Fraction(0), []
    for x in itertools.product(range(q), repeat=n):
        if len(set(x)) <= k:
            continue
        cnt = [0] * (len(links) + 1)
        for mask in range(1 << len(links)):
            got = {(s, r) for li, (s, r) in enumerate(links) if mask >> li & 1}
            decs = []
            for i in procs:
                vw = (x[i],) + tuple(x[j] if (j, i) in got else None for j in procs if j != i)
                sv = sorted({u for u in vw if u is not None})
                decs.append(sv[0] if len(sv) == 1 else rule[(i, vw)])
                assert decs[-1] in sv
            if len(set(decs)) > k:
                cnt[bin(mask).count("1")] += 1
        e = sum(c * p ** j * (1 - p) ** (len(links) - j) for j, c in enumerate(cnt))
        rows.append((x, cnt))
        worst = max(worst, e)
    return worst, rows


def main():
    n, q, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    os.makedirs("out", exist_ok=True)
    for spp in sys.argv[4:]:
        p = Fraction(spp)
        st, tv, rule, inputs, geo, W = run(n, q, k, p)
        worst, rows = evaluate(n, q, k, rule, p)
        print(f"n={n} q={q} k={k} p={p}: {st} OPT={Fraction(tv, W)} sim={worst} match={Fraction(tv, W) == worst}")
        st1 = run(n, q, k, p, cap=tv - 1)[0] if tv > 0 else "SKIP(OPT=0)"
        st0 = run(n, q, k, p, cap=tv)[0]
        print(f"   control t<=OPT-1: {st1} (must be INFEASIBLE); t<=OPT: {st0} (must be FEASIBLE or OPTIMAL)")
        for x, c in rows:
            print(f"   input {''.join(map(str, x))} fail counts {c}")
        with open(f"out/anonrule_n{n}_q{q}_k{k}_p{p.numerator}_{p.denominator}.json", "w", encoding="utf-8") as f:
            json.dump([[i, list(vw), v] for (i, vw), v in sorted(rule.items(), key=str)], f)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
