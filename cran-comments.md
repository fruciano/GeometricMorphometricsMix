## Test environments
* Local: R 4.5.1 on Windows 11 Pro
* GitHub Actions: ubuntu-latest (release, devel, oldrel), macOS-latest, windows-latest
* R-hub: R-devel/Fedora Linux 38, R-devel/Ubuntu 22.04.5 LTS


## R CMD check results

0 errors | 0 warnings | 1 note

## Notes
* This is a resubmission. Changes made in response to reviewer feedback:
  - Expanded the Description field to provide more detail about package functionality and implemented methods.
  - Added foundational references for the main methodological approaches in the Description field.
  - Replaced `\dontrun{}` with `\donttest{}` in examples for functions where examples are executable but take >5 seconds or require Suggests packages:
    - `Kmultparallel()` in R/Kmultparallel.R (requires phytools, MASS, mvMORPH, future; runs simulations)
    - `print.parallel_Kmult()` in R/Kmultparallel.R (depends on Kmultparallel output)
    - `plot.parallel_Kmult()` in R/Kmultparallel.R (depends on Kmultparallel output)
    - `summary.parallel_Kmult()` in R/Kmultparallel.R (depends on Kmultparallel output)

### Additional information
During local checks with vignette building enabled (Quarto + Pandoc installed), a cosmetic message appears once:

```
ERROR: Unknown command "TMPDIR=...". Did you mean command "create"?
```

This originates from an internal `system2()` call (in upstream tooling, not in this package) that invokes `quarto -V` while passing the temporary directory via the `env` argument on Windows, which Quarto interprets as a positional argument. The vignette build proceeds successfully (Pandoc 3.6.3; Quarto 1.8.25) and the message has no impact on the rendered vignettes. Setting `TMPDIR` in the environment prior to running the check suppresses it. No action required for CRAN.

The single NOTE reported is:

```
checking for future file timestamps ... NOTE
	unable to verify current time
```

This is an occasional Windows timing artifact (system clock / filesystem resolution); continuous integration does not reproduce additional issues, and there are no generated files with future timestamps. All other checks pass cleanly.
