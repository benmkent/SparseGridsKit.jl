# Spectral Sparse Grids
The sparse grid approximations are polynomial approximations, with a basis of Lagrange interpolation polynomials.
The approximation can be written in any suitable polynomial basis, for example we could use the Chebyshev polynomials.

The change of basis is implemented exploiting the `ApproxFun.jl` package.
The underlying interpolation polynomials used in `SparseGridsKit.jl` are constructed using `Fun`s.
To change basis we use an `ApproxFun.jl` expansion into coefficients and form a tensor product.
To promote sparsity this is done using a truncated Kronecker product.
The polynomial coefficients are assembled with corresponding polynomial degrees and spaces in a callable `SpectralSparseGridApproximation`.

## Examples
An adaptive sparse grid approximation for a Genz test function is created. 
```@example genz
using SparseGridsKit
n = 8
C = 1.0
W = 0.0
T = "exponentialdecay"
N = "gaussianpeak"
f = genz(n::Int, C::Float64, W::Float64, T::String, N::String)
domain = fill([-1,1],n)
# Approximate
(sg, f_on_Z) = adaptive_sparsegrid(f, domain, n)
f_sg = SparseGridApproximation(sg,f_on_Z)
```
The approximation `f_sg` is a representation formed as a linear combination of Lagrange interpolation polynomials.
This is converted to a spectral representation using [`convert_to_spectral_approximation`](@ref).
```@example genz
f_spectral = convert_to_spectral_approximation(sg, f_on_Z)
```
This still has an `n` dimensional domain.
```@example genz
f_spectral.dims
```
The representation of the tensor of polynomials is compactly stored as a `SparseVector`.
To index into this we require the storage dimension for each domain dimension.
```@example genz
f_spectral.expansiondimensions
```
This demonstrates that more polynomial terms are used in the parameter dimensions closer to one.
Intuitively, this matches the anisotropic nature of the test function (see [`genz`](@ref)).
The polynomial coefficients and the corresponding polynomial degrees are stored in sparse vectors.
```@example genz
f_spectral.coefficients, f_spectral.polydegrees
```
The spectral representation is the same function as the sparse grid representation.
```@example genz
x_test = [rand(n) for i in 1:100]
y_test = f_sg.(x_test)
y_spectral = f_spectral.(x_test)
all(isapprox(y_test, y_spectral; atol=1e-8))
```

The spectral sparse grid approximation also supports addition and subtraction.
```@example genz
f_spectral_2 = f_spectral + f_spectral
y_spectral_2 = f_spectral.(x_test)

all(isapprox(y_spectral_2, 2*y_spectral; atol=1e-8))
```
```@example genz
f_spectral_0 = f_spectral - f_spectral
y_spectral_0 = f_spectral.(x_test)

all(isapprox(y_spectral_0, zeros(size(y_spectral_0)); atol=1e-8))
```

## Function Reference
```@autodocs
Modules = [SparseGridsKit]
Pages   = ["spectralsparsegrids.jl"]
```
