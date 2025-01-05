using UMBridge
using SparseGridsKit
using Interpolations
using Plots

url = "http://stokes.imati.cnr.it:32716/openfoam"

# Set up a model by connecting to URL and selecting the "forward" model
model = UMBridge.HTTPModel("forwardcp", url)

n=1
k=5
mi = create_smolyak_miset(n,k)
sg = create_sparsegrid(mi)
#sg = create_sparsegrid(mi, rule=linear, knots=lejapoints)

Z = get_grid_points(sg)
# Z = [0.00000000000000e+000
# 1.00000000000000e+000
# -1.00000000000000e+000
# -577.350268922479e-003
# 577.350268922479e-003
# -839.943480719486e-003
# 839.943480719486e-003
# -276.062583588279e-003
# 276.062583588279e-003
# -943.018892243908e-003
# 943.018892243908e-003
# -434.632833901420e-003
# 434.632833901420e-003
# -732.417690451028e-003
# 732.417690451028e-003]
# Z = sort(Z)

# Model evaluation with configuration parameters
config = Dict("Fidelity" => 4, "res_tol" => 1e-10);
configlofi = Dict("Fidelity" => 1, "res_tol" => 1e-4);

linmap(x,a,b) = a + x*(b-a)
#map(z) = [linmap(0.5*(z+1), 2.34,23.4), 34.6]
map(z) = [23.4, linmap(0.5*(z+1), 3.46,34.6)]

cp = [UMBridge.evaluate(model, [map(z[1])], config) for z in Z]
cplofi = [UMBridge.evaluate(model, [map(z[1])], configlofi) for z in Z]

nhifi=660
nlofi=350

interp_at(x,cp;n=nhifi) = interpolate((cp[1][1:n],), cp[2][1:n],Gridded(Linear()))(x)

x = 0.8
cp_at_x = [interp_at(x,c) for c in cp] 
cplofi_at_x = [interp_at(x,c;n=nlofi) for c in cplofi]

zvec = [z[1] for z in Z]
zplot = range(-1, 1,500)
cp_poly = interpolate_on_sparsegrid(sg, cp_at_x, zplot)
cplofi_poly = interpolate_on_sparsegrid(sg, cplofi_at_x, zplot)
zvectoplot = linmap.(0.5*(zvec.+1), 3.46,34.6)
zplottoplot = linmap.(0.5*(zplot.+1), 3.46,34.6)
p = plot()
scatter!(p,zvectoplot, cp_at_x, marker=:x, ms=5, mc=:blue)
plot!(p, zplottoplot, cp_poly, label="hifi", color=:blue)
scatter!(p, zvectoplot, cplofi_at_x, marker=:x, ms=5, mc=:red)
plot!(p, zplottoplot, cplofi_poly, label="lofi", color=:red)
plot!(p,xlabel!("Flow speed U"))
plot!(p,ylabel!("C_p at x/c=0.8"))
plot!(p,title!("C_p for varied flow speed, fixed jet=23.4"))
display(p)