## Test environments
* Local: R 4.5.1 on Windows 11 Pro
* GitHub Actions: ubuntu-latest (release, devel, oldrel), macOS-latest, windows-latest


## R CMD check results

0 errors | 0 warnings | 2 notes


## Changes since the CRAN version (0.6.0.1)

This is an update to an existing CRAN package. Key changes:

* Made `Kmultparallel()` and its S3 methods (`print.parallel_Kmult()`, `plot.parallel_Kmult()`, `summary.parallel_Kmult()`) defunct. These were deprecated in the previous development version (0.6.1.0). All implementation code and examples have been removed; calling these functions now signals an informative error directing users to `geomorph::physignal()`.
* Improved the permutation loop in `pls()` to reduce its memory footprint.
* Added optional parallelism to `pls()` via the `future`/`future.apply` framework.
* Moved `ape` from `Depends` to `Suggests` since no remaining exported function uses `ape` directly.


## Notes

### NOTE 1: CRAN incoming feasibility — URL returning 403

```
Found the following (possibly) invalid URLs:
  URL: https://www.physalia-courses.org/courses-workshops/course22/
    From: README.md
    Status: 403
    Message: Forbidden
```

This URL in README.md is valid and accessible in a browser. The Physalia Courses server blocks automated/bot HTTP requests, producing a 403 status for `R CMD check`. The link is correct and intentional.


### NOTE 2: Future file timestamps

```
checking for future file timestamps ... NOTE
  unable to verify current time
```

This is an occasional network/system artifact on Windows when the check cannot reach a time server. There are no files with future timestamps in the package. Continuous integration on other platforms does not reproduce this NOTE.
