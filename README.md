# Optimal one-round set agreement for three processes through stochastic channels

Yitzchak Shmalo, Einstein Institute of Mathematics, The Hebrew University of Jerusalem.

Code, outputs and paper source for the paper of the same title. Three synchronous processes with identifiers
exchange one round of messages, and each of the six directed links delivers its message independently with
probability p. Each process holds one of Q >= 3 input values and decides a value it has seen; the execution fails
when the three decisions are distinct (2-set agreement). The paper determines the least worst-input failure
probability of a deterministic algorithm for every p and every Q >= 3. It has four polynomial pieces, attained by
four explicit rules C, S, L and m, with breakpoints 1 - 1/sqrt 2, the real root r = 0.352201... of
2x^3 - 4x^2 + 4x - 1, and 1/2.

All counts are exact integers and all reported values exact rationals.

## Layout

- `paper/main.tex`, `paper/main.pdf` - the paper.
- `code/` - the programs. Run them from the repository root.
- `out/` - the outputs the paper uses.

## Programs

- `code/env650.py` - the lower bound for Q = 3 (Proposition 8 of the paper): enumerates the 1728 source rules,
  counts failing delivery patterns for each, removes dominated count vectors, and computes the lower envelope of
  the remaining polynomials with exact real roots (SymPy). Also prints the envelope at the grid points used below.
- `code/leader_check.py` - direct simulation of the four rules C, S, L and m on every input with three distinct
  values for Q = 3, 4, 5, 6, compared with the counts in Table 1 of the paper. Standard library only.
- `code/setag3.py` - the check of Section 7: the exact optimum over all deterministic algorithms at a rational p,
  as a CP-SAT model with integer weights, with no reduction to source rules. At each p it also runs the controls
  t <= OPT - 1 (must be infeasible) and t <= OPT (must be feasible), and re-evaluates the returned algorithm by
  direct simulation in rational arithmetic. It writes the returned algorithm to `out/setrule_*.json` in the
  current directory (these files are not shipped).
- `code/milp3.py` - helper module imported by `setag3.py` (processes, links, views). Run on its own it is the
  mixed integer program for consensus from the companion paper; it is not used for any output here.

## What each output is and how it was made

| Output | Command | What it shows |
| --- | --- | --- |
| `out/env650.txt` | `python code/env650.py > out/env650.txt` | 1728 source rules, 201 distinct count vectors, 16 undominated; the lower envelope has four pieces with the counts of C, S, L and m, and the source rule attaining each |
| `out/leader_check.txt` | `python code/leader_check.py > out/leader_check.txt` | the counts of C, S, L and m on all 6, 24, 60 and 120 inputs with three distinct values for Q = 3..6; `ALL OK` |
| `out/out_grid.txt` | `/usr/bin/time -l python code/setag3.py 3 3 2 1/100 1/20 1/10 1/5 1/4 1/3 2/5 1/2 3/5 2/3 3/4 4/5 9/10 19/20 99/100 > out/out_grid.txt 2>&1` | CP-SAT optimum for Q = 3 at 15 values of p, with controls and the per-input counts of the returned algorithm; the last lines are the resource report of `time -l` (macOS) |
| `out/out_q4.txt` | `/usr/bin/time -l python code/setag3.py 3 4 2 2/5 9/20 1/5 1/3 3/5 > out/out_q4.txt 2>&1` | the same for Q = 4 (24 inputs with three distinct values) at five values of p |
| `out/envelope_vs_cpsat.txt` | the bash command below | the envelope of `env650.txt` against the CP-SAT optimum of `out_grid.txt` at the 15 points: p, envelope value, p, CP-SAT value; all `SAME` |

`out/envelope_vs_cpsat.txt` joins the two other outputs (bash, from the repository root):

```bash
paste -d' ' \
  <(grep '^  p=' out/env650.txt | sed -E 's/^  p=([^ ]+) OPT=([^ ]+)$/\1 \2/') \
  <(grep '^n=3 q=3' out/out_grid.txt | sed -E 's/^.* p=([^:]+): [A-Z]+ OPT=([^ ]+) .*$/\1 \2/') \
  | awk '{print $0, ($1==$3 && $2==$4) ? "SAME" : "DIFF"}' > out/envelope_vs_cpsat.txt
```

The arguments of `setag3.py` are the number of processes, the number of input values Q, k, and the values of p.

The optimum, the controls and the simulated value are determined by the model. The per-input counts printed under
each p are those of the algorithm the solver happens to return, and where several algorithms are optimal they can
change between runs: a rerun of the Q = 3 grid gave the same file except at p = 1/2, where L and m tie and the
solver returned a different optimal algorithm (same optimum, different counts on non-worst inputs).

## Requirements

Python 3.9.6 on an Apple M6 Mac mini, with OR-Tools 9.15.6755, SciPy 1.13.1, NumPy 2.0.2 and SymPy 1.14.0.
`setag3.py` needs OR-Tools and SciPy, `env650.py` needs SymPy, and `leader_check.py` needs only the standard
library. `env650.py` runs in about 8 seconds; each grid of `setag3.py` took under 15 seconds and under 600 MB.

## References

The model is due to P. Fraigniaud, B. Patt-Shamir and S. Rajsbaum, Coordination through stochastic channels,
DISC 2025, LIPIcs 356, 32:1-32:19, doi:10.4230/LIPIcs.DISC.2025.32. The broadcast variant is in their DISC 2026
paper, Consensus with stochastic broadcast, LIPIcs 397, 18:1-18:17, doi:10.4230/LIPIcs.DISC.2026.18.

The companion paper on consensus for three processes in the same model, with its code, is at
https://github.com/yspennstate/one-round-consensus-three-processes.
