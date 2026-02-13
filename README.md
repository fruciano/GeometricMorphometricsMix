# GeometricMorphometricsMix

Miscellaneous geometric morphometric functions

This package provides a diverse set of functions for geometric morphometric analyses, including shape integration, disparity, association, and more. It is designed for researchers, students, and practitioners in geometric morphometrics.

Originally developed as a resource for students in the [Geometric Morphometrics](https://www.physalia-courses.org/courses-workshops/course22/) course (Physalia Courses), these functions have been complemented over the years with many new functions, and are now distributed as an R package for broader accessibility and ease of use.

## Installation
From Github, using devtools

```
library(devtools)
install_github("fruciano/GeometricMorphometricsMix", build_vignettes=TRUE)
```

## Citation
Currently, there is no formal publication describing the package (but this may be available in the future). Please, use citation() to identify the most appropriate current citation of *GeometricMorphometricsMix*.
Importantly, each function provides reasonable comments and references to the methods, which can (and should) be used when using the functions in the package (in addition to the usual citation to the package itself).


## Documentation and Help

Each function in the package is documented with usage examples and references.


## Currently available functions
- **adjRand_test** : Permutation test of the adjusted Rand index, which quantifies the agreement between two partitions (e.g., classification of specimens obtained using different methods)
- **BTailTest** : Comparison of disparity between two groups
- **critical_angle** : Critical angle for the test of difference between two vectors
- **disparity_resample** : Obtain resampling estimates (bootstrap or rarefaction) of disparity statistics
- **disparity_test** : Comparison of disparity between two groups
- **dist_mean_boot** : Bootstrap estimates of the distance between group means
- **EscoufierRV** : Computation of Escoufier RV, which quantifies levels of association between blocks of variables (Escoufier 1973 - Biometrics)
- **Kmultparallel** : Parallelised computation of Adams' Kmult (useful for distributions of trees, see Fruciano et al. 2017 - Ecology and Evolution) - deprecated
- **LM_relativepos_check** : Check relative position of landmarks, comparing it to a reference specimen (useful to identify switched landmarks and similar in raw data)
- **parallel_analysis** : Perform parallel analysis (useful to choose a number of principal components for dimensionality reduction)
- **pls** : Perform partial least squares (PLS) analysis
- **pls_major_axis** : For a partial least squares (PLS) analysis, compute the major axis of PLS scores, project data onto the axis and compute predicted shapes for the major axis scores thus obtained.
- **ProjectOrthogonal** : Projects data to subspace orthogonal to a given vector (for use in allometric correction, see Burnaby 1966 - Biometrics; Rohlf & Bookstein 1987 - Systematic Zoology; for use in dealing with measurement error, see Valentin et al. 2008 - Journal of Fish Biology; Fruciano 2016 - Development Genes and Evolution)
- **rescale_by_landmark_distance** : Convenience function to rescale configurations of landmarks based on a vector of inter-landmark distances
- **reversePCA** : Simple function to obtain the original variables (e.g., shape) from PC scores (and mean)
- **rotate_landmarks**: Apply a user-defined rotation of a landmark configuration about the origin
- **repeated_measures_test** : Test of difference between two repeated measures
- **RVrarefied** : Computation of rarefied estimates of Escoufier RV (see Fruciano et al. 2013 - Plos One)
- **RVcomparison** : Permutation test for the difference in Escoufier RV (see Fruciano et al. 2013 - Plos One)
- **scaled_variance_of_eigenvalues** : Compute estimates of the scaled variance of eigenvalues (a commonly used measure of integration)
- **TestOfAngle** : Test of angle between two vectors, optionally allowing "flipping" of one of the two

# Function Groups by Analysis Type

| Analysis Type | Functions | Description |
|--------------|-----------|-------------|
| **Analysis of Disparity** | `BTailTest` | Traditional comparison of disparity between two groups |
| | `disparity_resample` | Resampling estimates of disparity |
| | `disparity_test` | Comparison of disparity between two groups |
| **Integration/Association/Modularity** | `EscoufierRV` | Computation of Escoufier RV coefficient, quantifying association between blocks of variables |
| | `RVrarefied` | Computation of rarefied estimates of Escoufier RV to account for sample size |
| | `RVcomparison` | Permutation test for the difference in Escoufier RV between datasets |
| | `scaled_variance_of_eigenvalues` | Compute estimates of the scaled variance of eigenvalues (measure of integration) |
| | `pls` | Perform partial least squares (PLS) analysis for covariation between blocks |
| | `pls_major_axis` | Compute major axis of PLS scores and predicted shapes along this axis |
| **Vector/Angle Analysis** | `critical_angle` | Critical angle for the test of difference between two vectors |
| | `TestOfAngle` | Test of angle between two vectors, with optional "flipping" of one vector |
| **Group Comparison and Classification** | `adjRand_test` | Permutation test of the adjusted Rand index, quantifying agreement between two partitions |
| | `dist_mean_boot` | Bootstrap estimates of the distance between group means |
| | `repeated_measures_test` | Test of difference between two repeated measures |
| **Data Transformation and Manipulation** | `ProjectOrthogonal` | Projects data to subspace orthogonal to a given vector (for allometric correction) |
| | `rescale_by_landmark_distance` | Rescale configurations of landmarks based on inter-landmark distances |
| | `reversePCA` | Obtain original variables from PC scores and mean |
| | `rotate_landmarks` | Apply a user-defined rotation of a landmark configuration |
| **Phylogenetic Comparative Analysis** | `Kmultparallel` | Parallelised computation of Adams' Kmult (for distributions of trees) - deprecated |
| **Quality Control and Diagnostics** | `LM_relativepos_check` | Check relative position of landmarks against a reference specimen |
| **Dimensionality Reduction** | `parallel_analysis` | Perform parallel analysis to determine number of principal components to retain |


## Badges
<!-- badges: start -->
  [![R-CMD-check](https://github.com/fruciano/GeometricMorphometricsMix/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fruciano/GeometricMorphometricsMix/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

