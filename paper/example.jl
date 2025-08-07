# Example for use in paper.md

using SparseGridsKit, Plots, LaTeXStrings
n = 2
C = 1.0
W = 0.0
T = "quadraticdecay"
N = "gaussianpeak"
f = genz(n, C, W, T, N)

# Adaptive approximation
(sg, f_on_Z) = adaptive_sparsegrid(f, n)

# Plot MI set and Sparse Grid
miset = get_mi_set(sg)
x     = getindex.(get_mi(miset), 1)
y     = getindex.(get_mi(miset), 2)
plot(
  scatter(x,y, xlabel=L"\beta_1", ylabel=L"\beta_2", legend=:none,
            title="Multi-Index Set", aspect_ratio=:equal),
  plot(sg; legend=:none, xlabel=L"y_1", ylabel=L"y_2", aspect_ratio=:equal)
)

# Plot approximation
plot(
    SparseGridApproximation(sg,f_on_Z);
    seriestype=:surface,
    title="",
    xlabel=L"y_1",
    ylabel=L"y_2",
    zlabel=L"u(\vec{y})"
)

# Demonstrate interpolation and quadrature
f_0 = interpolate_on_sparsegrid(sg, f_on_Z, [[0,0]])
integral = integrate_on_sparsegrid(sg,f_on_Z)
println("f([0,0]) = $f_0, Integral(f) = $integral")

# Demonstrate structures
f_sga = SparseGridApproximation(sg,f_on_Z)
f_spectral = convert_to_spectral_approximation(sg, f_on_Z)
f_sga_0 = f_sga([0,0])
f_spectral_0 = f_spectral([0,0])
println("f_sga([0,0]) = $f_sga_0, f_spectral(f) = $f_spectral_0")
