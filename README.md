# Verifier for “A proof of Borcea's quadratic variance conjecture”

This repository accompanies the paper *A proof of Borcea's quadratic variance
conjecture* by Ludovick Bouthat and Shaun Fallat (2026). It checks the finite
part of the proof of Theorem 1.1, for every degree 4 ≤ n ≤ 2047, by the method
of Sections 5–6 of the paper. It uses only the five restrictions of
Proposition 5.1:

- (1) 0 ≤ v ≤ s < 1
- (2) n s (1 + A) ≥ n^{2/m} m
- (3) |√A − √r| ≤ L
- (4) |Im μ| ≤ L
- (5) max{A + r − L², 1 + Ar − kv} ≤ b ≤ 2√(Ar)

Here m = n − 1, k = (m − 1)/m and L = min{(1 − s)√(mA), √((1 − s)(A − 1/m))}.

## Run

Python ≥ 3.10, standard library only; nothing to install. From this folder:

```text
python verify_all.py --workers 2
```

(On Windows, `py verify_all.py --workers 2` also works.) The run takes a few
minutes with 2 workers; more workers make it faster, and the worker count does
not affect the result. It prints one line per degree or degree block and ends
with this summary (only `seconds` varies):

```text
{
  "status": "PASS",
  "degrees": [
    4,
    2047
  ],
  "restrictions_used": [
    "(1)",
    "(2)",
    "(3)",
    "(4)",
    "(5)"
  ],
  "joint": {
    "sha256": "8fb676db95c9e8d513d3fe2bc3ecc165444cee13e38c876abd834f169ab15c41",
    "regions": 150845
  },
  "blocks": {
    "sha256": "69fe29be74825b9dba78e08c637e3ff57603cc5224f54dec17453929e59a76a8",
    "regions": 27429
  },
  "seconds": 150.0
}
```

The exit code is 0 exactly when the verification passes. A per-degree report is
written to `reports/report.json`; the copy shipped here is the output of the
reference run.

The script runs four steps:

1. Self-tests:
   - interval arithmetic encloses exact results;
   - quadrature weights are nonnegative and sum to the exact moments;
   - the coefficient identities behind Lemma 5.6 and Proposition 5.8 hold exactly;
   - the quartic-lemma constants check;
   - the restriction code is sound on random feasible points;
   - forged, unresolved or incomplete certificates are rejected;
   - a box around a known obstruction is not certified.
2. Replay of every box of `certificates/joint.json` (three parameters, degrees 4–26).
3. Replay of every region of `certificates/blocks.json` (two parameters, degree blocks 27–2047).
4. A check that every degree 4–2047 is covered exactly once.

Every box is rebuilt from the root and its bisections. Every leaf is recomputed.
Leaf labels are never trusted.

## Checking the files

`MANIFEST.sha256` lists the SHA-256 digest of every file except the report. On
Linux, `sha256sum -c MANIFEST.sha256`; on macOS, `shasum -a 256 -c MANIFEST.sha256`;
on any system with Python:

```text
python -c "import hashlib; [print('OK  ' if hashlib.sha256(open(f,'rb').read()).hexdigest()==h else 'FAIL', f) for h,f in (l.split() for l in open('MANIFEST.sha256'))]"
```

The digests of the two certificates are also printed in Appendix A of the paper.

## Files (what the proof relies on)

| File | Content | Paper |
|---|---|---|
| `borcea/interval.py` | integer interval arithmetic, outward rounding, 2^96 scale | §6.2 |
| `borcea/quadrature.py` | chord rule and the two rational meshes | Lemma 6.2, §6.3 |
| `borcea/restrictions.py` | Prop. 5.1 (1)–(5); K from Lemma 2.3 | Prop. 5.1, Lemma 2.3 |
| `borcea/endpoint.py` | shifted endpoint test and endpoint duality | Prop. 5.8, Cor. 5.9 |
| `borcea/joint.py` | three-parameter exclusions, 4 ≤ n ≤ 26 | §5, §6.1 |
| `borcea/blocks.py` | two-parameter exclusions for degree blocks | §6.4 |
| `borcea/certificate.py` | tree replay and coverage | Lemma 6.1 |
| `verify_all.py`, `tests/selftest.py` | driver and self-tests | |

`tools/refine.py` generated the certificates. It is **not** part of the proof: its output
is only accepted after a full replay by `verify_all.py`.

## How the restrictions are used

All of this is in `restrictions.py`.

- **(2)** is used in three ways:
  - `violates_2` excludes a region. For a block, the base is bounded by
    n₋ s(1+A)/(n₋−1) and raised to the power n₊−1 (§6.4).
  - In a single degree, s ≥ c_n/(1+A) narrows s.
  - The joint checker re-tests (2) after A has been narrowed.
- **(3)** is used to:
  - narrow A to [(√r − L)₊², (√r + L)²] (three-parameter checker only);
  - bound √r by [a − L, a + L].
- **(5)** is used as follows:
  - Its second lower bound against b ≤ 2a√r is the quadratic (A+k)x² − 2ax + 1 − ks ≤ 0
    in x = √r. So √r lies between its roots x∓ = (a ∓ √Δ)/(A+k), with Δ = k((A+k)s − 1),
    and Δ < 0 excludes the region.
  - Its first lower bound against the upper bound is equivalent to (3).
  - Lower bounds for b come from both lower bounds of (5), evaluated at the smallest
    admissible r: A + r − L², As + r + (1−s)/m, 1 + Ar − k(s − r), 2A(1 − L/a) (at
    r = (a−L)², by (3)), and 2a x₋.
  - A region is excluded when the lower bound for b exceeds 2√(Ar) (or 2aU in the
    two-parameter checker, where U is the largest admissible |μ|).
- **(1)** gives √r ≤ √s; the coordinates keep 0 ≤ v ≤ s ≤ 1. The strict s < 1 is what
  makes the centered boundary criteria strict (K√s < 1).
- **(4)** is used only in the dual test, to cap |Im μ| ≤ L in the enclosure of μ
  (the code encloses aμ, with |Im(aμ)| ≤ aL).

Other facts used, all proved elsewhere in the paper:

- A ≤ m (Lemma 2.1: the search space);
- A ≥ 1/m, a consequence of Lemma 2.2(4), which justifies writing (A − 1/m)₊ in L;
- |D(1)|² ≤ kv (Lemma 2.2(3));
- Π|p(0)| ≤ K = √A (s(n−A)/m)^{m/2} ≤ s^{m/2} (Lemma 2.3);
- β ≤ 1 and the convexity of β^ρ (Lemma 3.1, Lemma 6.2);
- the direct, centered, coefficient and endpoint estimates of Sections 3 and 5;
- the quartic local lemma with radius 1/20, recognized by 2401 v ≤ s (Proposition 5.10, §5.6).

## Exclusion names in the reports

`derivative-product` (2), `restriction-1`, `restriction-3`, `restriction-5`, `domain` (A > m),
`negative-endpoint` / `negative-interior` / `negative-D` (β ≥ 0, |D(1)|² ≥ 0), `direct`,
`centered-coefficient`, `centered-coefficient-boundary`, `centered-integral`,
`centered-integral-boundary`, `centered-monotone` (6.1), `shifted:T`, `dual:T`, `quartic-local`.

## Certificate format

- **`borcea-joint-v2`:** one tree per degree 4–26. The root is
  [0, n−1] × [0, 1] × [0, 1] in (A, e, f), with e = √(1−s) (η in the paper) and f = √(v/s).
  Internal nodes are `{"axis": 0|1|2, "children": [lo, hi]}`. Leaves are
  `{"test": name, "parts": q}` with 8 ≤ q ≤ 512, using the uniform mesh.
- **`borcea-blocks-v2`:** consecutive blocks covering 27–2047. Each block has a root
  [0, n₊−1] × [0, 1] in (A, e) and splits `"A"`, `"e"` or `"n"` (the degree interval).
  Leaves are `{"test", "parts"}` with 8 ≤ q ≤ 4096. A block may set
  `"quadrature": "geometric"`.

## License

- Code (everything except `certificates/` and `reports/`): MIT License, see `LICENSE`.
- Certificates and report: CC BY 4.0, see `LICENSE-DATA.md`.

## Citation

If you use this verifier or its certificates, please cite the paper and this
archive; `CITATION.cff` contains the metadata (GitHub shows it under
“Cite this repository”).
