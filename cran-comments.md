## Test environments

* Local: R 4.5.1 on Windows 11 Pro
* GitHub Actions: ubuntu-latest (release, devel, oldrel), macOS-latest, windows-latest

## R CMD check results

0 errors | 0 warnings | 1 note

## Changes since the CRAN version (0.6.1.1)

This is an update to an existing CRAN package, submitted as version 0.6.1.3. Key changes:

* Added median output to `disparity_resample()` and made the median the default point estimate shown by its plotting method.
* Expanded unit tests for existing functions and added regression tests for the fixes below.
* Fixed conversion of the second three-dimensional input array in `repeated_measures_test()`.
* Fixed missing bootstrap confidence-interval endpoints in `scaled_variance_of_eigenvalues()` when the number of bootstrap replicates is small, by bounding the rounded interval indices to valid rows.
* Added `withr` to Suggests to restore temporary state in tests.

## Notes

### NOTE: Future file timestamps

```
checking for future file timestamps ... NOTE
  unable to verify current time
```

The local Windows check could not verify the current time, producing this NOTE.
