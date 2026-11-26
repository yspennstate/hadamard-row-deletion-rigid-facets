"""Theorem 1.1(a) and Remark 3.3 on normalized Paley type I Hadamard matrices of orders 12 and 20. For every choice of
second deleted row, no word of the truncated cube outside the columns is within the diameter of all of them
(exhaustive). The three-row pattern count (N + s1 s2 s3 T)/8 is compared with direct counts. Needs numpy."""
import numpy as np, itertools

def paley(q):
    sq = {(x * x) % q for x in range(1, q)}
    chi = lambda a: 0 if a % q == 0 else (1 if a % q in sq else -1)
    Q = np.array([[chi(j - i) for j in range(q)] for i in range(q)])
    S = np.zeros((q + 1, q + 1), int)
    S[0, 1:] = 1; S[1:, 0] = -1; S[1:, 1:] = Q
    H = S + np.eye(q + 1, dtype=int)
    return H

def normalize(H):
    H = H * H[:, :1]          # first column +1
    H = H * H[:1, :]          # first row +1
    return H

for q in (11, 19):
    H = normalize(paley(q)); N = len(H)
    assert (H @ H.T == N * np.eye(N)).all() and (H[0] == 1).all() and (H[:, 0] == 1).all()
    r = N // 2
    L = N - 2
    words = np.array(list(itertools.product((1, -1), repeat=L)), dtype=np.int8) if L <= 18 else None
    bad = 0
    for s in range(1, N):
        keep = [i for i in range(N) if i not in (0, s)]
        C = H[keep].T.astype(np.int8)                        # N columns of length L
        d = (L - C.astype(int) @ C.T.astype(int)) // 2
        diam = d.max()
        assert diam == r, (q, s, diam)
        colset = {tuple(c) for c in C}
        ok = np.ones(len(words), bool)
        for c in C:                                            # distance <= r  <=>  dot >= L - 2r
            ok &= (words.astype(np.int16) @ c.astype(np.int16)) >= L - 2 * r
        extra = [w for w in words[ok] if tuple(w) not in colset]
        bad += bool(extra)
    print(f"Paley order {N}: two-deletion choices with an addable word: {bad} of {N-1} (diameter r={r} each)")
    # three deletions: pattern count formula (N + s1 s2 s3 sum a b t)/8 against direct counts
    mism = 0
    for a, b, t in itertools.permutations(range(1, N), 3):
        tri = (H[a] * H[b] * H[t]).sum()
        for s1, s2, s3 in itertools.product((1, -1), repeat=3):
            cnt = int(((H[a] == s1) & (H[b] == s2) & (H[t] == s3)).sum())
            if 8 * cnt != N + s1 * s2 * s3 * tri:
                mism += 1
    print(f"   pattern-count formula mismatches over all ordered triples: {mism}")
