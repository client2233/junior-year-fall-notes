#!/usr/bin/env python3
NEG = float("-inf")


def pr(F, x, y, title):
    print(f"\n-- {title} --")
    print("      " + "".join(f"{c:>6}" for c in ("-" + y)))
    for i, row in enumerate(F):
        lbl = "-" if i == 0 else x[i - 1]
        print(f"{lbl:>4} " + "".join("   -inf" if v == NEG else f"{v:>6}" for v in row))


def gotoh(x, y, match=1, mismatch=-1, d=-3, e=-1):
    n, m = len(x), len(y)
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    Ix = [[NEG] * (m + 1) for _ in range(n + 1)]
    Iy = [[NEG] * (m + 1) for _ in range(n + 1)]
    Pm = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    Px = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    Py = [[[] for _ in range(m + 1)] for _ in range(n + 1)]
    M[0][0] = 0
    for i in range(1, n + 1):
        Ix[i][0] = d + (i - 1) * e
        Px[i][0] = ["Ix"] if i > 1 else ["M"]
    for j in range(1, m + 1):
        Iy[0][j] = d + (j - 1) * e
        Py[0][j] = ["Iy"] if j > 1 else ["M"]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if x[i - 1] == y[j - 1] else mismatch
            cand = [(M[i - 1][j - 1], "M"), (Ix[i - 1][j - 1], "Ix"), (Iy[i - 1][j - 1], "Iy")]
            best = max(v for v, _ in cand)
            Pm[i][j] = [k for v, k in cand if v == best]
            M[i][j] = best + s
            cand = [(M[i - 1][j] + d, "M"), (Ix[i - 1][j] + e, "Ix")]
            best = max(v for v, _ in cand)
            Px[i][j] = [k for v, k in cand if v == best]
            Ix[i][j] = best
            cand = [(M[i][j - 1] + d, "M"), (Iy[i][j - 1] + e, "Iy")]
            best = max(v for v, _ in cand)
            Py[i][j] = [k for v, k in cand if v == best]
            Iy[i][j] = best
    return M, Ix, Iy, Pm, Px, Py


x, y = "AAT", "ACACT"
M, Ix, Iy, Pm, Px, Py = gotoh(x, y)
pr(M, x, y, "M  (x[i] ~ y[j])")
pr(Ix, x, y, "Ix (x[i] ~ gap)")
pr(Iy, x, y, "Iy (y[j] ~ gap)")
n, m = len(x), len(y)
best = max(M[n][m], Ix[n][m], Iy[n][m])
print(f"\nbest = {best}  M={M[n][m]} Ix={Ix[n][m]} Iy={Iy[n][m]}")

# enumerate all optimal alignments
out = []


def walk(state, i, j, ax, ay):
    if state == "M" and i == 0 and j == 0:
        out.append(("".join(reversed(ax)), "".join(reversed(ay))))
        return
    if state == "M":
        for k in Pm[i][j]:
            s = match if x[i - 1] == y[j - 1] else mismatch
            prev = {"M": M, "Ix": Ix, "Iy": Iy}[k][i - 1][j - 1]
            if prev + s != M[i][j]:
                continue
            walk(k, i - 1, j - 1, ax + [x[i - 1]], ay + [y[j - 1]])
    elif state == "Ix":
        if i == 0:
            return
        for k in Px[i][j]:
            base = M[i - 1][j] if k == "M" else Ix[i - 1][j]
            if base + (d_ := (-3 if k == "M" else -1)) != Ix[i][j]:
                continue
            if i == 1 and j == 0 and k == "M":
                out.append(("".join(reversed(ax)), "".join(reversed(ay))))
                continue
            walk(k, i - 1, j, ax + [x[i - 1]], ay + ["-"])
    else:
        if j == 0:
            return
        for k in Py[i][j]:
            base = M[i][j - 1] if k == "M" else Iy[i][j - 1]
            if base + (d_ := (-3 if k == "M" else -1)) != Iy[i][j]:
                continue
            if j == 1 and i == 0 and k == "M":
                out.append(("".join(reversed(ax)), "".join(reversed(ay))))
                continue
            walk(k, i, j - 1, ax + ["-"], ay + [y[j - 1]])


for st, val in (("M", M[n][m]), ("Ix", Ix[n][m]), ("Iy", Iy[n][m])):
    if val == best:
        walk(st, n, m, [], [])

print("\noptimal alignments:")
for a, b in sorted(set(out)):
    print(f"x: {a}\ny: {b}\n")
