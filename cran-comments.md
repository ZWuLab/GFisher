# cran-comments for GFisher 0.3.1

## Test environment
- R 4.3.3 on x86_64-conda-linux-gnu
- check command: `R CMD check --as-cran`
- LaTeX present: full manual check ran.

## R CMD check results
- 0 ERROR | 0 WARNING | 3 NOTE

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

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable
```
