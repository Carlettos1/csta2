# Library design document

[design.tex](design.tex) is the complete English source;
[design.pdf](design.pdf) is the compiled document. It covers all six crates,
physical conventions, algorithms, code examples, supported ranges, numerical
contracts, performance costs, and limitations of the current local implementation.

Build from the workspace root:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=docs/build docs/design.tex
cp docs/build/design.pdf docs/design.pdf
```

The document uses standard TeX Live packages: `babel`, `lmodern`, `microtype`,
`geometry`, AMS math, `booktabs`, `longtable`, `tabularx`, `xcolor`, `listings`,
`enumitem`, `fancyhdr`, PGF/TikZ, `hyperref`, and `bookmark`. It requires neither
shell escape nor external images or bibliography processing. TeX Live 2026 and
`latexmk` were available when the document was created.

Validate the complete embedded Rust examples against the local library:

```sh
python3 docs/check_examples.py --offline
python3 docs/check_examples.py --offline --run
```

Omit `--offline` if dependencies need downloading. The checker requires Rust
1.98.0 or newer, checks default and checkpoint configurations, and gives every
example a finite runtime timeout. It creates temporary source files and uses
`target/design-examples` for Cargo artifacts. It does not modify workspace members
or the workspace lockfile. Trait, manifest, and shell excerpts are not standalone
programs and are excluded from extraction.

The examples validate usage and selected invariants. Their small stochastic
budgets do not establish convergence or publication-quality physical estimates.

Validation on 6 September 2026: the 51-page PDF built without TeX warnings,
representative pages were visually inspected, and all 22 complete Rust examples
compiled and ran successfully. The source-path references were also checked.
