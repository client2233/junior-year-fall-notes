#!/usr/bin/env python3
"""NW / SW / Gotoh worked examples for the notes."""

NEG = float("-inf")


def grid(F, x, y, title):
    print(f"\n-- {title} --")
    print("      " + "".join(f"{c:>5}" for c in ("-" + y)))
    for i, row in enumerate(F):
        lbl = "-" if i == 0 else x[i - 1]
        cells = []
        for v in row:
            cells.append("  -inf" if v == NEG else f"{v:>5}")
        print(f"{lbl:>4} " + "".join(cells))


def nw_all(x, y, match=2, mismatch=-1, gap=-2):
    n, m = len(x), len(y)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap
    P = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if x[i - 1] == y[j - 1] else mismatch
            c = [(F[i - 1][j - 1] + s, "D"), (F[i - 1][j] + gap, "U"), (F[i][j - 1] + gap, "L")]
            best = max(v for v, _ in c)
            P[i][j] = [d for v, d in c if v == best]
            F[i][j] = best
    paths = []

    def walk(i, j, ax, ay):
        if i == 0 and j == 0:
            paths.append(("".join(reversed(ax)), "".join(reversed(ay))))
            return
        if i == 0:
            walk(0, j - 1, ax + ["-"], ay + [y[j - 1]])
            return
        if j == 0:
            walk(i - 1, 0, ax + [x[i - 1]], ay + ["-"])
            return
        for d in P[i][j]:
            if d == "D":
                walk(i - 1, j - 1, ax + [x[i - 1]], ay + [y[j - 1]])
            elif d == "U":
                walk(i - 1, j, ax + [x[i - 1]], ay + ["-"])
            else:
                walk(i, j - 1, ax + ["-"], ay + [y[j - 1]])

    walk(n, m, [], [])
    uniq = sorted(set(paths))
    return F, uniq, F[n][m]


def sw(x, y, match=1, mismatch=-1, gap=-2):
    n, m = len(x), len(y)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    bi = bj = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if x[i - 1] == y[j - 1] else mismatch
            c = [(F[i - 1][j - 1] + s, "D"), (F[i - 1][j] + gap, "U"), (F[i][j - 1] + gap, "L"), (0, "0")]
            best = max(v for v, _ in c)
            P[i][j] = [d for v, d in c if v == best]
            F[i][j] = best
            if best > F[bi][bj]:
                bi, bj = i, j
    paths = []

    def walk(i, j, ax, ay):
        if F[i][j] == 0:
            paths.append(("".join(reversed(ax)), "".join(reversed(ay))))
            return
        for d in P[i][j]:
            if d == "0":
                paths.append(("".join(reversed(ax)), "".join(reversed(ay))))
            elif d == "D":
                walk(i - 1, j - 1, ax + [x[i - 1]], ay + [y[j - 1]])
            elif d == "U":
                walk(i - 1, j, ax + [x[i - 1]], ay + ["-"])
            else:
                walk(i, j - 1, ax + ["-"], ay + [y[j - 1]])

    walk(bi, bj, [], [])
    return F, sorted(set(paths)), F[bi][bj], (bi, bj)


print("#" * 70)
print("# 1) NW homework: x=GAATTCAGTTA  y=GGATCGA   match2/mismatch-1/gap-2")
F, paths, sc = nw_all("GAATTCAGTTA", "GGATCGA")
grid(F, "GAATTCAGTTA", "GGATCGA", "F matrix")
print(f"\nscore = {sc}")
for a, b in paths:
    print(f"x: {a}\ny: {b}\n")

print("#" * 70)
print("# 2) NW class ex: x=ACGCATCA y=ACTGATTCA  match2/mismatch-1/gap-2")
F, paths, sc = nw_all("ACGCATCA", "ACTGATTCA")
grid(F, "ACGCATCA", "ACTGATTCA", "F matrix")
print(f"\nscore = {sc}")
for a, b in paths:
    print(f"x: {a}\ny: {b}\n")

print("#" * 70)
print("# 3) NW class AAAC/AGC  match1/mismatch-1/gap-2 (slide s=1/-1/-2)")
F, paths, sc = nw_all("AAAC", "AGC", match=1, mismatch=-1, gap=-2)
grid(F, "AAAC", "AGC", "F matrix")
print(f"\nscore = {sc}")
for a, b in paths:
    print(f"x: {a}\ny: {b}\n")

print("#" * 70)
print("# 4) SW demo: x=GAATTCAGTTA y=GGATCGA  match1/mismatch-1/gap-2")
F, paths, sc, pos = sw("GAATTCAGTTA", "GGATCGA")
grid(F, "GAATTCAGTTA", "GGATCGA", "F matrix (local)")
print(f"\nmax score = {sc} at F{pos}")
for a, b in paths:
    print(f"x: {a}\ny: {b}\n")
