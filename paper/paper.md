---
title: 'SparseGridsKit.jl: Adaptive single- and multi-fidelity sparse grid approximation in Julia'
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
    # equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
    corresponding: true # (This is how to denote the corresponding author)
affiliations:
 - name: CNR-IMATI, Pavia, Italy
   index: 1
date: 20 May 2025
bibliography: paper.bib
---

# Summary
Approximation of functions with high dimensional domains is important for modern scientific and engineering problems.
An example of this is constructing surrogate models for quantities of interest for high dimensional parametrised PDE problems.
These surrogate models are constructed to give computationally cheap yet accurate approximations that can be used in applications such as uncertainty quantification, optimisation and parameter estimation [@UQHandbook2017].
Surrogates may be constructed with global polynomial approximation on the parameter space and a common approach is the use of *sparse grid* approximation techniques.
In particular, sparse grid polynomial interpolation techniques allow a practitioner to approximate solutions to parametric problems in a non-intrusive way using existing numerical solvers.

`SparseGridsKit.jl` provides a Julia toolbox to manually and adaptively construct sparse grid polynomial approximations [@julia].
Interpolation and quadrature routines allow evaluation and integration of the surrogate models.
Multi-fidelity approximation via the multi-index stochastic collocation algorithm is also possible [@HajiAli2016] [@Jakeman2019] [@Piazzola2022].
Approximations can be represented either in a basis of Lagrange interpolation polynomials or in a basis of spectral-type polynomials.

# Statement of need
Sparse grid approximation is a well developed methodology and is featured in many survey articles and textbook chapters, e.g. [@Bungartz2004],[@LeMaitre2010],[@Schwab2011],[@Cohen2015],[@Sullivan2015].
The need for sparse grid surrogate modelling is demonstrated by its use in many applications, from simpler elliptic and parabolic PDEs to complex practical engineering problems e.g.\ [@Piazzola2021],[@Piazzola2022],[@Li2024].
The `SparseGridsKit.jl` implementation offers a rich set of features to enable this.

Specifically, `SparseGridsKit.jl` is a Julia implementation of sparse grid approximation methods.
This offers

- native Julia implementation of adaptive sparse grid approximation functionality,
- dynamical typing, allowing surrogate models to map input parameters to any Julia type offering vector space operations.

Existing sparse grid approximation packages in Julia include [`Tasmanian.jl`](https://github.com/floswald/Tasmanian.jl), wrapping the [Tasmanian library](https://github.com/ORNL/Tasmanian), [`AdaptiveSparseGrids.jl`](https://github.com/jacobadenbaum/AdaptiveSparseGrids.jl) and [`DistributedSparseGrids.jl`](https://github.com/baxmittens/DistributedSparseGrids.jl).
`SparseGridsKit.jl` offers a more complete set of functionality, with close resemblance to the popular `Sparse Grids MATLAB Kit` [@Piazzola2024].

Other popular software packages implementing sparse grid approximation include:

- `Sparse Grids MATLAB Kit`: A MATLAB package on which the `SparseGridsKit.jl` is loosely based [@Piazzola2024].
- `spinterp`: A MATLAB toolbox for sparse grid interpolation [@spinterp] (no longer maintained).
- `UQLab`: A broad MATLAB uncertainty quantification toolkit [@Marelli2014].
- `PyApprox`: A Python package for high-dimensional approximation [@PyApprox].
- `Dakota`: A C++ library for optimisation and surrogate modelling [@Dakota].
- `UQTk`: A collection of C++/Python uncertianty quantification tools including sparse grid quadrature [@Debusschere2015].
- `Tasmanian`,`SG++`,: C++ sparse grid approximation implementations with wrappers for many popular software languages [@stoyanov2015tasmanian] [@pflueger10spatially].

`SparseGridsKit.jl` offers specific Julia toolkit with minimal complexity for fast algorithm development and prototyping.

# `SparseGridsKit.jl` Features
The main features are outlined below:

- One dimensional knots and quadrature rules.
- Multi-index set construction and manipulation.
- Combination technique sparse grid approximations including evaluation and quadrature routines.
- Adaptive sparse grid approximation construction based on the ubiquitous Gerstner-Griebel dimensional adaptive algorithm [@Gerstner2003]. 
- Adaptive multi-fidelity approximation via the Multi-Index Stochastic Collocation (MISC) algorithm [@HajiAli2016] [@Jakeman2019] [@Piazzola2022].
- Conversion to and from Polynomial Chaos / spectral polynomial series representation.
- Limited support for surrogate model differentiation via automatic differentiation.

The functionality described above is tested and documented with examples included in the repository.

# Acknowledgements
The author has been supported by the project 202222PACR “Numerical approximation of uncertainty quantification problems for PDEs by multi-fidelity methods (UQ-FLY)", funded by European Union -- NextGenerationEU.

# References