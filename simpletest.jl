using SparseGridsKit
# using GLMakie
using DistributedSparseGrids
using StaticArrays 

# GLMakie.activate!()
nparams = 2
fun1(x,ID="") =  (1.0-exp(-1.0*(abs(2.0 - (x[1]-1.0)^2.0 - (x[2]-1.0)^2.0) +0.01)))/(abs(2-(x[1]-1.0)^2.0-(x[2]-1.0)^2.0)+0.01)

function sparse_grid(N::Int,pointprops; nlevel=6,RT=Float64,CT=Float64)
	# define collocation point
	CPType = CollocationPoint{N,CT}
	# define hierarchical collocation point
	HCPType = HierarchicalCollocationPoint{N,CPType,RT}
	# init grid
	asg = init(AHSG{N,HCPType},pointprops)
	#set of all collocation points
	cpts = Set{HierarchicalCollocationPoint{N,CPType,RT}}(collect(asg))
	# fully refine grid nlevel-1 times
	for i = 1:nlevel-1
		union!(cpts,generate_next_level!(asg))
	end
	return asg
end

(sg, f_on_Z) = adaptive_sparsegrid(fun1, nparams, maxpts=10_000, proftol=1e-3)
pcl = precompute_lagrange_integrals(7)
integral_result = integrate_on_sparsegrid(sg, f_on_Z, pcl)

asg = sparse_grid(2, @SVector [1,1])
init_weights!(asg, fun1)
for i = 1:20
	cpts = generate_next_level!(asg, 1e-3, 20)
	init_weights!(asg, collect(cpts), fun1)
end

nptsperdim = 200
xs = LinRange(-1.0, 1.0, nptsperdim)
zvals = [[x,y] for x in xs, y in xs];
@time zs1 = interpolate_on_sparsegrid(sg, f_on_Z, zvals[:]);
@time zs2 = [interpolate(asg,[x,y]) for x in xs, y in xs];

# cpts_x_1 = [DistributedSparseGrids.coord(hcpt,1) for hcpt in asg]
# cpts_x_2 = [DistributedSparseGrids.coord(hcpt,2) for hcpt in asg]
# cpts_x_3 = zeros(length(cpts_x_1))


# fig = Figure();
# ax3d1 = Axis3(fig[1,1]);
# ax3d2 = Axis3(fig[1,2]);
# surface!(ax3d1, xs, xs, reshape(zs1, nptsperdim, nptsperdim));
# cpts_sg = get_grid_points(sg)
# scatter!(ax3d1, map(x->x[1], cpts_sg), map(x->x[2], cpts_sg), zeros(length(cpts_sg)), markersize=2)

# surface!(ax3d2, xs, xs, zs2);
# scatter!(ax3d2, cpts_x_1, cpts_x_2, cpts_x_3, markersize=2)

# fig