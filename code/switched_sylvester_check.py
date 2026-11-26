"""Theorem 1.1(b): the switched Sylvester matrix H(x,y) = (-1)^(x.y + x1 x2 [y in E]), E = e3 + span(e1,e2), of order
N = 2^m. Delete rows 0, e1, e2 and k-3 further rows other than t = e1+e2. Check that the truncated columns have
diameter N/2 and that the four columns in E with their t-entry flipped are new words, within N/2 of every column and of
each other. The same test on the plain Sylvester matrix (no switch) is run as a control and must fail.
Columns are bit-packed integers over the kept rows. Run: python3 switched_sylvester_check.py (a few minutes)."""
import sys


def run(m, k, switched=True):
    N = 1 << m
    bit = lambda v, i: (v >> i) & 1
    inE = lambda y: bit(y, 2) == 1 and (y >> 3) == 0
    entry = lambda x, y: (bin(x & y).count('1') + (bit(x, 0) * bit(x, 1) * inE(y) if switched else 0)) & 1  # 1 is -1
    t = 0b11
    deleted = [0, 1, 2] + [x for x in range(3, N) if x != t][: k - 3]
    if m <= 8:  # Hadamard and normalized, checked exhaustively
        rows = [[entry(x, y) for y in range(N)] for x in range(N)]
        assert all(v == 0 for v in rows[0]) and all(r[0] == 0 for r in rows)
        for i in range(N):
            for j in range(i + 1, N):
                assert sum(1 if a == b else -1 for a, b in zip(rows[i], rows[j])) == 0
    keep = [x for x in range(N) if x not in deleted]
    pos_t = keep.index(t)
    cols = []
    for y in range(N):
        w = 0
        for p, x in enumerate(keep):
            if entry(x, y):
                w |= 1 << p
        cols.append(w)
    r = N // 2
    pc = lambda v: bin(v).count('1')
    colset = set(cols)
    diam = max(pc(a ^ b) for a in cols for b in cols)
    new = [cols[y] ^ (1 << pos_t) for y in range(N) if inE(y)]
    ok_new = all(w not in colset for w in new) and len(set(new)) == 4
    ok_all = all(max(pc(w ^ c) for c in cols) <= r for w in new)
    ok_mut = max(pc(a ^ b) for a in new for b in new) <= r
    return diam, r, ok_new and ok_all and ok_mut


for m in range(3, 11):
    for k in (3, 4, 5):
        if (1 << m) < (1 << k):
            continue
        d, r, ok = run(m, k)
        _, _, ctrl = run(m, k, switched=False)
        print(f"N={1 << m:5d} k={k} diameter={d} (N/2={r})  four words adjoinable: switched={ok}  plain Sylvester={ctrl}")
        assert d == r and ok and not ctrl
        sys.stdout.flush()
print("all checks passed")
