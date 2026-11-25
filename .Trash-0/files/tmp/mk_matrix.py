#!/usr/bin/env python3
import sys

TMP = "/home/LittleWatter/Documents/Obsidian Vault/大三上/.workbuddy/tmp"


def load(name):
    rows = []
    with open(f"{TMP}/{name}.txt", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(line.split())
    header = rows[0]
    order = header
    mat = {}
    for r in rows[1:]:
        mat[r[0]] = dict(zip(order, [int(v) for v in r[1:]]))
    return order, mat


def md(name, order, mat):
    out = ["| " + " | ".join([""] + order) + " |", "|" + "---|" * (len(order) + 1)]
    for a in order:
        out.append("| **" + a + "** | " + " | ".join(str(mat[a][b]) for b in order) + " |")
    return "\n".join(out)


for n in sys.argv[1:]:
    order, mat = load(n)
    print(f"#### {n}")
    print(md(n, order, mat))
    print()
