#!/usr/bin/env python3
"""Needleman-Wunsch trace output for the class homework example."""

MATCH, MISMATCH, GAP = 2, -1, -2


def nw(x, y, match=MATCH, mismatch=MISMATCH, gap=GAP):
    n, m = len(x), len(y)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap
    P = [[None] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if x[i - 1] == y[j - 1] else mismatch
            cands = [(F[i - 1][j - 1] + s, "D"), (F[i - 1][j] + gap, "U"), (F[i][j - 1] + gap, "L")]
            best = max(c[0] for c in cands)
            P[i][j] = [c[1] for c in cands if c[0] == best]
            F[i][j] = best
    # traceback (take first preference D > U > L)
    i, j = n, m
    ax, ay, path = [], [], []
    while i > 0 or j > 0:
        path.append((i, j))
        if i == 0:
            ax.append("-"); ay.append(y[j - 1]); j -= 1
            continue
        if j == 0:
            ax.append(x[i - 1]); ay.append("-"); i -= 1
            continue
        if i > 0 and j > 0 and "D" in P[i][j]:
            ax.append(x[i - 1]); ay.append(y[j - 1]); i -= 1; j -= 1
        elif i > 0 and "U" in P[i][j]:
            ax.append(x[i - 1]); ay.append("-"); i -= 1
        else:
            ax.append("-"); ay.append(y[j - 1]); j -= 1
    return F, P, "".join(reversed(ax)), "".join(reversed(ay)), F[n][m]


def show(x, y):
    print(f"x = {x}   (n={len(x)})")
    print(f"y = {y}   (m={len(y)})")
    print(f"match={MATCH} mismatch={MISMATCH} gap={GAP}\n")
    F, P, ax, ay, sc = nw(x, y)
    m = len(y)
    header = "      " + "".join(f"{'-' if c == '-' else c:>5}" for c in ("-" + y))
    print(header)
    for i, row in enumerate(F):
        label = "-" if i == 0 else x[i - 1]
        print(f"{label:>4} " + "".join(f"{v:>5}" for v in row))
    print(f"\nscore F(n,m) = {sc}")
    print(f"x: {ax}")
    print(f"y: {ay}")
    ident = sum(1 for a, b in zip(ax, ay) if a == b and a != "-")
    print(f"identical pairs = {ident}, alignment length = {len(ax)}")


print("=" * 60)
print("HW (deck 03-1 slide 75):  GAATTCAGTTA  vs  GGATCGA")
print("=" * 60)
show("GAATTCAGTTA", "GGATCGA")

print()
print("=" * 60)
print("In-class example (deck 03-1 ss 55-64):  AAAC vs AGC")
print("=" * 60)
show("AAAC", "AGC")

print()
print("=" * 60)
print("In-class example (deck 03-1 ss 67-73): x=ACTGATTCA  y=ACGCATCA")
print("=" * 60)
show("ACTGATTCA", "ACGCATCA")
