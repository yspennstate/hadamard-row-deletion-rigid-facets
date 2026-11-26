# Row deletion in normalized Hadamard matrices and rigid facets of Vietoris-Rips complexes of hypercubes

Yitzchak Shmalo, Einstein Institute of Mathematics, The Hebrew University of Jerusalem.

Briggs, Feng and Wells (Facets in the Vietoris-Rips complexes of hypercubes, Exp. Math. 35 (2026) 388-396, doi:10.1080/10586458.2025.2474034, arXiv:2408.01288) conjectured that for every fixed k the columns of a normalized Hadamard matrix of large order, with its first k rows deleted, form a maximal set of their diameter in the Hamming cube, and proved this for k at most 1. We prove it for k = 2 and every order N >= 4. For every k >= 3 the conjecture fails: at each order 2^m >= 2^k a switched Sylvester matrix gives a counterexample, and four new words can be adjoined to it without increasing the diameter. We also answer their question whether the number f'(n;r) of dimensions in which the Vietoris-Rips complex VR(Q_n;r) has a rigid facet grows exponentially in r. It is at least 2^(r-6)+1 for r >= 6 and n >= r+2, and a theorem of Bukh and Saatashvili (arXiv:2507.10828) bounds it by (4-10^(-10))^r for large r. Finally, the smallest facet of VR(Q_n;r) and the smallest rigid facet have the same size once n >= 2^r.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/switched_sylvester_check.py`: Theorem 1.1(b). Builds the switched Sylvester matrix of order N = 8 to 1024, checks that it is a normalized Hadamard matrix (exhaustively up to N = 256), deletes k = 3, 4, 5 rows and checks the four adjoinable words; the plain Sylvester matrix is run as a control and fails at every (N, k). Run `python3 switched_sylvester_check.py`; a few minutes.
- `code/paley_check.py`: Theorem 1.1(a) and Remark 3.3 on the normalized Paley matrices of orders 12 and 20, by exhaustive search over the truncated cube. Needs numpy; about a minute.
- `code/fibre_check.py`: Lemma 4.1 (the sets A and B) and the antipodal-fibre construction of Theorem 1.2 for L = 1, 2, 3. Under a minute.
- `SHA256SUMS`: sha256 of each file above.
