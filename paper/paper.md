---
title: 'SparseGridsKit.jl: A Julia sparse grids approximation implementation'
tags:
  - Julia
  - numerical approximation
  - sparse grids
  - uncertainty quantification
authors:
  - name: Benjamin M. Kent
    orcid: 0000-0003-4968-7993
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
    corresponding: true # (This is how to denote the corresponding author)
affiliations:
 - name: CNR-IMATI, Pavia, Italy
   index: 1
date: 20 January 2025
bibliography: paper.bib
---

# Summary
Approximation of functions with high dimensional domains is an important tool for modern scientific and engineering modelling.
Global polynomial approximation is often justifiable and a common approach is the use of *sparse grid* approximation techniques.
In particular, sparse grid polynomial interpolation techniques allow a practitioner to use existing numerical solvers and approximate solutions to parametric formulations of such problems in a non-intrusive way.
In the context of uncertainty quantification, postprocessing of such an approximation allows a user to easily and cheaply extract expected values, variances, Sobel indices and may other useful descriptive properties. 

# Statement of need

`SparseGridsKit.jl` offers a Julia language implementation of sparse grid approximation methods.
A user can construct grids offering refinement in different domain dimensions by hand crafting the underlying *multi-index set* structure or using the Gerstner--Griebel [@Gerstner2001] style adaptive construction `adaptive_sparsegrid`.
The package is designed with a focus on uncertainty quantification.
For a grid and corresponding set of function evaluations a user can interpolate the approximation to points in the function domain, compute weighted integrals of the approximation over thh function domain, and compute weighted integrals of the square of the approximation over the function domain.

This package is strongly related to the `Sparse Grids MATLAB Kit` [@Piazzola2024].
It is hoped that `SparseGridsKit.jl` can exploit the improved efficiency of a Julia implementation, and also exploit the code flexibility encouraged by Julia.
For example, for interpolation and integration the function evaluations must simply return objects that satisfy vector space operations.
An example is given in the documentation using `ApproxFun` objects.
Computation of weighted $L^2$ type integrals is also possible if the user provides an implmentation of a norm for the vector space objects.

# Related Software

# Sparse Grid Approximation

# Benchmarks

# Acknowledgements
