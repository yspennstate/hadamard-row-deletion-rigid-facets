"""Lemma 4.1 and the construction in the proof of Theorem 1.2: the seven-bit sets A and B (integer encoding of Section 4)
are rigid and maximally diameter-5 in Q_7 with cross distances at most 6, and the antipodal-fibre sets over Q_L,
L = 1, 2, 3, are rigid and maximally diameter-(L+5) for every choice pattern, with the sizes 27*2^L + 2j."""
import itertools
A = [6,16,25,39,41,43,56,62,63,74,75,76,78,83,91,93,95,99,101,104,106,112,114,119,122,124,126]
B = [16,24,37,38,41,42,45,46,55,56,62,68,75,77,78,92,95,98,100,105,106,108,112,115,120,124,125,126]
pc = lambda v: bin(v).count('1')

def props(S, n):
    Sset = set(S)
    diam = max(pc(a ^ b) for a in S for b in S)
    rigid = all(max(pc(a ^ b) for b in S) == diam for a in S)
    maximal = all(max(pc(x ^ b) for b in S) > diam for x in range(1 << n) if x not in Sset)
    return diam, rigid, maximal

print('A', len(A), props(A, 7)); print('B', len(B), props(B, 7))
print('cross diameter', max(pc(a ^ b) for a in A for b in B))
for L in (1, 2, 3):
    pairs = [(p, p ^ ((1 << L) - 1)) for p in range(1 << L) if p < (p ^ ((1 << L) - 1))]
    sizes = set()
    for choice in itertools.product((A, B), repeat=len(pairs)):
        S = [(p << 7) | w for (p, q), F in zip(pairs, choice) for base in (p, q) for w in F for p in [base]]
        d, rg, mx = props(S, L + 7)
        assert d == L + 5 and rg and mx, (L, d, rg, mx)
        sizes.add(len(S))
    print(f'L={L}: all {2**len(pairs)} choices maximal, rigid, diameter {L+5}; sizes {sorted(sizes)}')
