---
title: 'SparseGridsKit.jl: Adaptive sparse grid approximation in Julia'
tags:
  - numerical approximation
  - surrogate modelling
  - sparse grids
  - multi-fidelity
  - uncertainty quantification
  - Julia
authors:
  - name: Benjamin M. Kent
    orcid: 0000-0003-4968-7993
    equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
    corresponding: true # (This is how to denote the corresponding author)
affiliations:
 - name: CNR-IMATI, Pavia, Italy
   index: 1
date: 20 January 2025
bibliography: paper.bib
---

# Summary
Approximation of functions with high dimensional domains is an important task in modern scientific and engineering modelling.
Surrogate models are often constructed to give computationally cheap yet accurate approximations that can be used in applications such as uncertainty quantification, optimisation, parameter estimation.
Surrogates may use global polynomial approximation and a common approach is the use of *sparse grid* approximation techniques.
In particular, sparse grid polynomial interpolation techniques allow a practitioner to approximate solutions to parametric formulations of such problems in a non-intrusive way using existing numerical solvers.

`SparseGridsKit.jl` provides tools to manually and adaptively construct sparse grid polynomial approximations.
Interpolation and quadrature routines allow evaluation and integration of the surrogate models.
Multi-fidelity approximation via the multi-index stochastic collocation algorithm is also possible.
Approximations can be represented either in a basis of Lagrange interpolation polynomials or in a basis of spectral-type polynomials.

# Statement of need
Surrogate modelling is an important theme in many science and engineering applications.
In particular, sparse grid approximation is a well developed methodology and is featured in many survey articles and textbook chapters, e.g. [@Bungartz2004,LeMaitre2010,@Schwab2011,@Cohen2015,Sullivan2015].
The need for sparse grid surrogate modelling is demonstrated by its use in many applications, from simpler elliptic and parabolic PDEs to complex practical engineering problems e.g.\ [@Piazzola2021,@Li2024,@Piazzola2022].
The `SparseGridsKit.jl` implementation offers a rich set of features to enable this.

Specifically, a Julia implementation of sparse grid approximation methods is offered by `SparseGridsKit.jl`.
This offers
  - native algorithm implementation in Julia of adaptive sparse grid approximation functionality,
  - dynamical typing, allowing surrogate models mapping parameters $\vec{y}$ to any type offering vector space operations.
Existing sparse grid approximation packages in Julia include [`Tasmanian.jl`](https://github.com/floswald/Tasmanian.jl), wrapping the [Tasmanian library](https://github.com/ORNL/Tasmanian), [`AdaptiveSparseGrids.jl`](https://github.com/jacobadenbaum/AdaptiveSparseGrids.jl) and [`DistributedSparseGrids.jl`](https://github.com/baxmittens/DistributedSparseGrids.jl).
`SparseGridsKit.jl` offers a more complete set of functionality, with close resemblance to the popular `Sparse Grids MATLAB Kit` [@Piazzola2024].

Other software packages implementing sparse grid approximation include:
  - `Sparse Grids MATLAB Kit` A MATLAB package on which the `SparseGridsKit.jl` is loosely based [@Piazzola2024],
  - `spinterp` A MATLAB toolbox for sparse grid interpolation [@spinterp],
  - `Dakota` A C++ library for optimisation and surrogate modelling [@Dakota].
  - `PyApprox` A Python package for high-dimensional approximation [@PyApprox],
`SparseGridsKit.jl` offers specific toolkit with minimal complexity for fast algorithm development in Julia.

# `SparseGridsKit.jl`
## Knots, Quadrature and Interpolation
Interpolation and quadrature methods are built out of one-dimensional knot functions which return knots $x$ and corresponding quadrature weights $w$.
Sparse grid constructions require sequences of approximation operators.
In this interpolation setting, we consider a sequence of interpolation operators $I^{\alpha}$ indexed by $\alpha\in\mathbb{N}$ which use sequentially growing sets of interpolation points.
This is achieved using the abstract `Points` and `Level` types in `SparseGridsKit.jl`.

Sparse grid methods extend one dimensional rules to many dimensions.
This originates in the work of Smolyak [@Smolyak1963] and developed for interpolation, quadrature and parametric PDEs [@Gerstner1998],[@Novak1999],[@Barthelmann2000],[@Babuska2004],[@Xiu2005],[@Nobile2008a].
In particular, we use the combination technique formulation in which the approximation is a cleverly chosen linear combination of tensor product interpolation operators.
<!-- The sparse grid interpolation operator is
$$
  I = \sum_{\underline{\alpha}\in A} c_{\alpha} \bigotimes_{i=1}^n I^{\alpha_i}
$$
where $c_{\underline{\alpha}}$ is the combination technique coefficient and $A$ is a set of multi-indices defining the approximation.. -->
For further details we refer to the vast literature, e.g. [@Piazzola2021].
Construction of an approximation using the sparse grid interpolation operator requires evaluation of the target function at a set of collocation points $Z$ which is implicitly defined.
Nesting of the underlying one-dimensional interpolation rules means that many grid points are coincident hence the *sparse* nature of the grid.

## Surrogate Models
For a sparse grid approximation with corresponding function evaluations $\{f(z)\}_{z\in Z}$, we can apply interpolation to approximate the function values at non-interpolation points.
This is implemented by the `interpolate_on_sparsegrid` function.
A `SparseGridApproximation` structure wraps the sparse grid approximation operator and a set of function evaluations which can be treated as a function to evaluate the approximation.
Similarly, the surrogate model can be integrated using the `integrate_on_sparsegrid` function.

## Adaptive Sparse Grids
A sparse grid can be constructed by using a user specified multi-index set.
Generally, the underlying function is unknown and we wish to adaptively construct the approximation space to capture the target function behaviour.
This is often achieved in a greedy iterative manner.
Adaptive sparse grid approximation is implemented as `adaptive_sparsegrid`.
This is based on the ubiquitous Gerstner-Griebel dimensional adaptive algorithm [@Gerstner2003].

## Other functionality
`SparseGridsKit.jl` includes functionality for multi-fidelity approximation via the multi-index stochastic collocation algorithm and limited support for differentiation via automatic differentiation.

The functionality described above is all tested and documented with examples.

# Simple Example


<!-- # Comparison of Features with Related Software
This package is strongly related to the `Sparse Grids MATLAB Kit` (SGMK) [@Piazzola2024].
The functionality of `SparseGridsKit.jl` is compared to the SGMK.
This is formally done in the automated testing.
This is not a complete comparison as we do not provide a complete reimplementation of the features in SGMK.

| Features                           | Common                                                                                                 | Differences                                                                                                                                                                                                                                                                                                             |
| ---------------------------------- | ------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Quadrature and Interpolation Rules | Equispaced<br>Gauss-Hermite<br>Gauss-Legendre<br>Weighted Leja<br>Clenshaw-Curtis                      | \`SparseGridsKit.jl\`<br>\- \`CustomPoints\` structure that allow user to supply a knots-weights function.<br>\- Leja points are computed online allowing custom weight functions.<br><br>\`SGMK\`<br>\- In built implementation of further knots-weights rules.                                                        |
| Adaptive Algorithm                 | Dimension-adaptive sparse grid algorithm.<br>Supports non-nested knots.                                | \`SparseGridsKit.jl\`<br>\- User input profit definitions.<br>\- Common single and multi-fidelity adaptive algorithm.<br>\- Support for any type implementation vector space operations.<br><br>\`SGMK\`<br>\- Provides multiple profit definitions<br>\- Implements a buffering strategy for high-dimensional problems |
| Parallelism                        |                                                                                                        | \`SparseGridsKit.jl\`<br>\- None built in<br><br>\`SGMK\`<br>\- Evaluations can use MATLAB Parallel Toolbox<br>\- Function evaluation recycling.                                                                                                                                                                        |
| Polynomial Chaos Expansion         | Convert polynomial approximations from Lagrange interpolation type basis to spectral polynomial basis. | \`SparseGridsKit.jl\`<br>\-<br><br>\`SGMK\`<br>\- Supports Legendre, Chebyshev, Hermite, Laguerre, Gen. Laguerre, Jacobi                                                                                                                                                                                                |
| Derivatives                        |                                                                                                        | \`SparseGridsKit.jl\`<br>\- No implementation of gradients.<br><br>\`SGMK\`<br>\- Finite difference computation of gradient and Hessian.<br>\- Global and local sensitivity via Sobol Indices.                                                                                                                          |
| Data Export and Visualization      | Plotting of sparse grids and sparse grid approximations.                                               | \`SparseGridsKit.jl\`<br>\- Plots are implemented as \`Recipies\` for use with \`Plots.jl\`.<br><br>\`SGMK\`<br>\- Export sparse grid points and corresponding weights to ASCII files.                                                                                                                                  | -->

# Acknowledgements

# References