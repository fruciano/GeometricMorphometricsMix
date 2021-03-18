# GeometricMorphometricsMix

Miscellaneous geometric morphometric functions

This repository contains a few useful and fairly diverse geometric morphometric functions which may be useful to others.
All the functions are reasonably well commented with information on usage and similar.
It started as a simple set of function that I distributed as help to the students of my course in [Geometric Morphometrics](https://www.physalia-courses.org/courses-workshops/course22/) (organized by Physalia Courses).
They are incorporated into an R package to facilitate distribution and use.

## Currently available functions
- **adjRand_test** : Permutation test of the adjusted Rand index, which quantifies the agreement between two partitions (e.g., classification of specimens obtained using different methods)
- **BTailTest** : Comparison of disparity between two groups
- **critical_angle** : Critical angle for the test of difference between two vectors
- **disparity_test** : Comparison of disparity between two groups
- **dist_mean_boot** : Bootstrap estimates of the distance between group means
- **EscoufierRV** : Computation of Escoufier RV, which quantifies levels of association between blocks of variables (Escoufier 1973 - Biometrics)
- **Kmultparallel** : Parallelised computation of Adams' Kmult (useful for distributions of trees, see Fruciano et al. 2017 - Ecology and Evolution)
- **LM_relativepos_check** : Check relative position of landmarks, comparing it to a reference specimen (useful to identify switched landmarks and similar in raw data)
- **parallel_analysis** : Perform parallel analysis (useful to choose a number of principal components for dimensionality reduction)
- **pls** : Perform partial least squares (PLS) analysis
- **pls_major_axis** : For a partial least squares (PLS) analysis, compute the major axis of PLS scores, project data onto the axis and compute predicted shapes for the major axis scores thus obtained.
- **ProjectOrthogonal** : Projects data to subspace orthogonal to a given vector (for use in allometric correction, see Burnaby 1966 - Biometrics; Rohlf & Bookstein 1987 - Systematic Zoology; for use in dealing with measurement error, see Valentin et al. 2008 - Journal of Fish Biology; Fruciano 2016 - Development Genes and Evolution)
- **rescale_by_landmark_distance** : Convenience function to rescale configurations of landmarks based on a vector of inter-landmark distances
- **reversePCA** : Simple function to obtain the original variables (e.g., shape) from PC scores (and mean)
- **rotate_landmarks**: Apply a user-defined rotation of a landmark configuration about the origin
- **rarefied_convex_hull** : Computation of rarefied estimates of n-dimensional convex hull volume (a measure of disparity/morphospace occupation)
- **rarefied_disparity** : Rarefied estimates of some measures of disparity/morphospace occupation
- **repeated_measures_test** : Test of difference between two repeated measures
- **RVrarefied** : Computation of rarefied estimates of Escoufier RV (see Fruciano et al. 2013 - Plos One)
- **RVcomparison** : Permutation test for the difference in Escoufier RV (see Fruciano et al. 2013 - Plos One)
- **scaled_variance_of_eigenvalues** : Compute estimates of the scaled variance of eigenvalues (a commonly used measure of integration)
- **TestOfAngle** : Test of angle between two vectors, optionally allowing "flipping" of one of the two


## Installation
From Github, using devtools

```
library(devtools)
install_github("fruciano/GeometricMorphometricsMix")
```

## Citation
Currently, and for the foreseeable future, there are no plans to publish a proper paper describing the package.
Each function, however, provides reasonable comments and references to the methods, which can (and should) be used when using the functions in the package (in addition to the usual citation to the package itself)

  <!-- badges: start -->
  [![R build status](https://github.com/fruciano/Geometric_morphometrics/workflows/R-CMD-check/badge.svg)](https://github.com/fruciano/Geometric_morphometrics/actions)
  <!-- badges: end -->
