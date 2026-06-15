# cran-comments for GFisher 0.3.0

## Test environment
- R 4.3.3 on x86_64-conda-linux-gnu
- check command: `R CMD check --as-cran --no-manual`
- NOTE: no LaTeX (pdflatex) was available in this environment, so the PDF manual was not built (`--no-manual`). The Rd files themselves pass all Rd checks; the PDF manual builds on a LaTeX-equipped machine (e.g. CRAN).

## R CMD check results
- 0 ERROR | 0 WARNING | 2 NOTE

- No ERRORs.

- No WARNINGs.

### NOTEs (verbatim from the check log)

```
* checking for future file timestamps ... NOTE
unable to verify current time
```

```
* checking compilation flags used ... NOTE
Compilation used the following non-portable flag(s):
  ‘-march=nocona’
```
