# Curtz's denominator conjecture: OEIS A189733, revision 12, COMMENTS, May 23, 2011

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang (National
University of Singapore; ChronoAI Pte Ltd, Singapore). Revision 2, November 3, 2026.
License: Creative Commons Attribution 4.0 International (CC BY 4.0),
https://creativecommons.org/licenses/by/4.0/.

## Result and sources

For every integer k >= 0, the six reduced denominators of the first row of
the rational array in https://oeis.org/A189733 are
(1, 3k+1, 2k+1, 3k+2, 1, 6k+6). Zero has denominator one.
The note follows Bing Gao's generating-function proof recorded in
https://oeis.org/history?seq=A189733, revision 30 (posted September 23, 2026;
contribution dated September 22). It also proves the numerator period, the
denominator recurrence for n >= 12, Colin Barker's empirical generating
function, and the inverse-binomial first-column identity. The original posing
edit is revision 12, COMMENTS (May 23, 2011); the conjecture moved to FORMULA
in revision 13. Additional corollaries give every array entry as a finite sum,
all numerator identities in https://oeis.org/A278718, and both Barker conjectures
and all three sequence identifications in https://oeis.org/A267068.
The second-kind companion `C(m,n)=2B(m,n+1)-B(m,n)` gives the denominator
formulas, generating function and reciprocal sum in https://oeis.org/A194767.
Its fourth-diagonal numerator formula needs `A146535(m+1)`, rather than
`A146535(m)`; at `m=2` the actual entry is `-5/12`.
For the unsigned-seed array in https://oeis.org/A189731 and
https://oeis.org/A242926, `u_n=(L_(n+1)-1)/(n+1)` with Lucas numbers `L_n`.
The odd denominator uses `A175386(r+1)`, not `A175386(r)`. A finite modular
certificate gives `L_705=1 (mod 705)`, so the integer at zero-based index 704
refutes the prime-only rank assertion. A general integer triangular
transformation proves the antidiagonal denominator-lcm formula in A242926
for both the unsigned and alternating-seed arrays. The source's informal
divisor and super-autosequence descriptions are not interpreted as precise identities.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, use a standard TeX distribution with amsmath, amsthm,
amssymb, lmodern, geometry, booktabs, array, longtable and hyperref:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is included in `note.tex`; no BibTeX database or figures are
required. A document using a separate bibliography may instead use the normal
pdflatex/bibtex/pdflatex/pdflatex sequence.

## Exact finite checks

Python 3.10 or later, using the standard library:

```sh
python3 verify.py --output verification.json
```

The checker computes the triangle from its boundary values and recurrence
with exact rational arithmetic. It compares first-row values, reduced
denominators, numerator periods, logarithmic coefficients, the denominator
generating function, and the stated recurrence through n=240. It checks the
full-square reflection and binomial identity through index 40, and the
polynomial moment identities through polynomial index 80. It checks the
numerator floor formula and order-five recurrence, the quotient tuple and
generating function through index 240, all three quotient identifications
through index 79, and 880 all-entry sums. Stored terms from the cited OEIS entries
are compared with their defining formulas, with the published offsets.
The second-kind denominator formulas and sequence identifications are checked
through index 240, and its first four diagonals through index 80. In total,
all 1,200 stored terms from eighteen entries are compared, with per-entry
counts in `verification.json`.
The unsigned row formula is compared with its recurrence through index 240;
the bisections through index 79 and antidiagonal sums through index 119 are
checked. Repeated matrix squaring and an independent Lucas recurrence both
check the residue at 705. The antidiagonal denominator-lcm invariant is checked
for both zero-diagonal arrays through antidiagonal 80.
All ranges and
counts appear in `verification.json`. These are finite checks; the note
contains the proofs for arbitrary indices.

## Lean statements

The pinned source is
https://github.com/the-omega-institute/trureturing/blob/20b69217fd7c02cbb8668d939dfa29ecc12a9716/D5/S1/Recurrence/Residue/OEISA189733DenominatorConjecture.lean.
In namespace `D5.S1.Recurrence.Residue.OEISA189733DenominatorConjecture`,
the public theorem `result` proves the universal denominator tuple and the
private theorem `triangle_six` proves the rational first-row tuple.
The generating-function proof, numerator-period conclusion, denominator
recurrence, denominator generating function, and full-square binomial
identity, all-entry finite sum, additional numerator and quotient identities,
and the second-kind formulas
are written proofs in the note. The formal definition `triangle m d` equals
`B(m,m+d)` and `a n` is its reduced denominator at `m=0`.

## Files

- `note.tex`: complete mathematical text and bibliography.
- `metadata.json`: title, creators, license and statement-source link.
- `verify.py`: exact-arithmetic checker, without third-party dependencies.
- `verification.json`: its deterministic output for the stated ranges.
- `source-terms.json`: cited sequence terms, offsets and version identifiers.
- `SHA256SUMS`: hashes of these scientific source files and this README.
