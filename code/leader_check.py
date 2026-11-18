"""Exact failure counts of four explicit rules on every input with three distinct values, Q = 3..6 (direct simulation).
C: lone message adopted; hearing both -> cyclic successor of own value among the three seen values.
S: lone message adopted; hearing both -> smallest heard value.
L: process 0 (leader) always keeps its value; others: lone message adopted; hearing both -> the leader's value.
m: smallest value seen.
Expect, for every input: C (1,0,3,5,3,3,1), S (1,0,3,5,4,1,0), L (1,2,2,2,1,0,0), m (1,3,3,1,0,0,0).
usage: python3 leader_check.py"""
import itertools

P = (0, 1, 2)
LINKS = [(s, r) for r in P for s in P if s != r]


def decide(rule, i, x, heard):
    own = x[i]
    vals = [x[j] for j in heard]
    if not vals:
        return own
    if rule == "m":
        return min([own] + vals)
    if rule == "L":
        if i == 0:
            return own
        if len(vals) == 1:
            return vals[0]
        return x[0]
    if len(vals) == 1:
        return vals[0]
    if rule == "S":
        return min(vals)
    seen = sorted([own] + vals)            # rule C: cyclic successor of own among the seen values
    return seen[(seen.index(own) + 1) % 3]


def counts(rule, x):
    c = [0] * 7
    for mask in range(64):
        got = {LINKS[t] for t in range(6) if mask >> t & 1}
        d = [decide(rule, i, x, [j for j in P if j != i and (j, i) in got]) for i in P]
        if len(set(d)) == 3:
            c[bin(mask).count("1")] += 1
    return tuple(c)


EXPECT = {"C": (1, 0, 3, 5, 3, 3, 1), "S": (1, 0, 3, 5, 4, 1, 0), "L": (1, 2, 2, 2, 1, 0, 0), "m": (1, 3, 3, 1, 0, 0, 0)}
bad = 0
for Q in range(3, 7):
    for rule, want in EXPECT.items():
        got = {counts(rule, x) for x in itertools.permutations(range(Q), 3)}
        ok = got == {want}
        bad += not ok
        print(f"Q={Q} rule {rule}: {len(list(itertools.permutations(range(Q), 3)))} inputs, count vectors {sorted(got)} {'OK' if ok else 'MISMATCH'}")
print("ALL OK" if not bad else f"{bad} MISMATCHES")
