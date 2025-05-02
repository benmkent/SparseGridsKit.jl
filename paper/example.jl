using SparseGridsKit, Plots

f(x) = x[1]^2 + x[1]*sin(x[2])

(sg, f) = adaptivesparsegrids()