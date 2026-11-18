"""Check that an anonymous rule file is a function of (own value, sorted heard values) only, and print its table.
usage: python3 anon_decode.py out/anonrule_n3_q4_k2_p2_5.json"""
import json
import sys

rule = json.load(open(sys.argv[1]))
table, clash = {}, 0
for i, vw, v in rule:
    key = (vw[0], tuple(sorted(h for h in vw[1:] if h is not None)))
    if key in table and table[key] != v:
        clash += 1
    table[key] = v
print(f"views {len(rule)}, distinct keys {len(table)}, clashes {clash} ({'ANONYMOUS' if clash == 0 else 'NOT anonymous'})")
lone = {k: v for k, v in table.items() if len(k[1]) == 1}
both = {k: v for k, v in table.items() if len(k[1]) == 2}
print("lone message (own, heard) -> keep/adopt:")
print("  " + "  ".join(f"{k[0]}{k[1][0]}:{'K' if v == k[0] else 'A'}" for k, v in sorted(lone.items())))
print("both heard (own; heard) -> decision:")
print("  " + "  ".join(f"{k[0]};{k[1][0]}{k[1][1]}:{v}" for k, v in sorted(both.items())))
