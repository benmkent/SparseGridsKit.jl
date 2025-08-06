@testset "Sparse Grids MATLAB Kit Tests" begin

using MAT, Downloads, FastGaussQuadrature

    function read_sgmk_mat(url)
        # Temporary file to store the download
        matfile = tempname() * ".mat"

        # Download the file
        Downloads.download(url, matfile)

        # Read the MAT file
        data = matread(matfile)
    return data
    end

    function compare_sg_to_sgmk(S,S_sgmk, Sr_sgmk)
        @test get_n_grid_points(S) == Sr_sgmk["size"]
        Z1 = sort(get_grid_points(S))
        Z2 = [z for z in eachcol(Sr_sgmk["knots"])]
        @test isapprox(Z1, Z2; atol=1e-8)
    end

    function compare_sg_to_sgmk(S,S_sgmk)
        Z_sgmk = [];
        for knots_matrix in S_sgmk["knots"]
            append!(Z_sgmk, eachcol(knots_matrix))
        end
        Z_sgmk = unique_with_tol(Z_sgmk)
        sort!(Z_sgmk)

        @test get_n_grid_points(S) == length(Z_sgmk)
        if get_n_grid_points(S) == length(Z_sgmk)
            Z = sort(get_grid_points(S))
            @test Z ≈ Z_sgmk
        end
    end

    function unique_with_tol(Z, tol=1e-12)
        unique_Z = []
        for z in Z
            if all(x -> !isapprox(x, z, atol=tol), unique_Z)
                push!(unique_Z, z)
            end
        end
        return unique_Z
    end

    @testset "Level to Knots" begin
    ii = [1, 2, 3, 4];
    m_lin = Linear().(ii);
    m_2step = TwoStep().(ii);
    m_doub = Doubling().(ii);
    m_trip = Tripling().(ii);
    # m_GK = lev2knots_GK.(ii);

    S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_lev2knots.mat")
    @test m_lin == vec(Integer.(S["m_lin"]))
    @test m_2step == vec(Integer.(S["m_2step"]))
    @test m_doub == vec(Integer.(S["m_doub"]))
    @test m_trip == vec(Integer.(S["m_trip"]))
    # @test m_GK == Integer(S["m_GK"])
    end

    @testset "Knots" begin
    n = 5; a = 1; b =4 ;
    # uniform pdf
    x_unif,w_unif         = GaussLegendrePoints([a,b])(n);
    x_CC,w_CC             = CCPoints([a,b])(n);
    x_leja,w_leja         = LejaPoints([a,b],false,:classic, z->sqrt(0.5))(n);
    # [x_sym_leja,w_sym_leja] = knots_leja(n,a,b,'sym_line');
    # [x_p_leja,w_p_leja]     = knots_leja(n,a,b,'p_disk');
    # [x_midp,w_midp]         = knots_midpoint(n,a,b);
    # [x_trap,w_trap]         = knots_trap(n,a,b);

    # normal pdf
    n = 9; mu = 0; sigma = 1;
    x_norm,w_norm                     = GaussHermitePoints()(n);
    # [x_GK,w_GK]                         = knots_GK(n,mu,sigma);
    # [x_norm_Leja,w_norm_Leja]           = knots_normal_leja(n,mu,sigma,'line');
    # [x_norm_sym_Leja,w_norm_sym_Leja]   = knots_normal_leja(n,mu,sigma,'sym_line');

    @debug "Skipping exponential PDF"
    # exponential pdf
    # n = 12; lambda = 1; 
    # [x_exp, w_exp]           = gausslaguerre(n,lambda);
    # # [x_exp_leja, w_exp_leja] = knots_exponential(n,lambda);

    @debug "Skipping gamma PDF"
    # # gamma pdf
    # n = 12; alpha = 1; beta = 2;
    # [x_gamma, w_gamma]           = knots_gamma(n,alpha,beta);
    # [x_gamma_leja, w_gamma_leja] = knots_gamma(n,alpha,beta);

    @debug "Skipping beta PDF"
    # beta pdf
    # n = 12; a = 1; b = 3; 
    # alpha = -0.5; beta = 0.5; 
    # [x_beta,w_beta]                   = knots_beta(n,alpha,beta,a,b);
    # [x_beta_leja,w_beta_leja]         = knots_beta_leja(n,alpha,beta,a,b,'line');
    # alpha = 1.5; beta = 1.5; % alpha = beta for symmetric points 
    # [x_beta_sym_leja,w_beta_sym_leja] = knots_beta_leja(n,alpha,beta,a,b,'sym_line');

    @debug "Skipping triangular PDF"
    # triangular pdf
    # n = 12; a = 0; b = 2;
    # [x_triang,w_triang] = knots_triangular_leja(n,a,b);

    S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_knots.mat")

    @test vec(x_unif) ≈ vec(S["x_unif"])
    @test vec(w_unif) ≈ vec(S["w_unif"])
    @test vec(x_CC) ≈ vec(S["x_CC"])
    @test vec(w_CC) ≈ vec(S["w_CC"])
    @test isapprox(vec(x_leja),vec(S["x_leja"]); rtol=1e-2) # These are known to be slightly different
    @test isapprox(vec(w_leja), vec(S["w_leja"]); rtol=1e-2)
    # @test x_sym_leja ≈ S["x_sym_leja"]
    # @test w_sym_leja ≈ S["w_sym_leja"]
    # @test x_p_leja ≈ S["x_p_leja"]
    # @test w_p_leja ≈ S["w_p_leja"]
    # @test x_midp ≈ S["x_midp"]
    # @test w_midp ≈ S["w_midp"]
    # @test x_trap ≈ S["x_trap"]
    # @test w_trap ≈ S["w_trap"]
    @test vec(x_norm) ≈ vec(S["x_norm"])
    @test vec(w_norm) ≈ vec(S["w_norm"])
    # @test x_GK ≈ S["x_GK"]
    # @test w_GK ≈ S["w_GK"]
    # @test x_norm_Leja ≈ S["x_norm_Leja"]
    # @test w_norm_Leja ≈ S["w_norm_Leja"]
    # @test x_norm_sym_Leja ≈ S["x_norm_sym_Leja"]
    # @test w_norm_sym_Leja ≈ S["w_norm_sym_Leja"]
    # @test x_exp ≈ S["x_exp"]
    # @test w_exp ≈ S["w_exp"]
    # @test x_exp_leja ≈ S["x_exp_leja"]
    # @test w_exp_leja ≈ S["w_exp_leja"]
    # @test x_gamma ≈ S["x_gamma"]
    # @test w_gamma ≈ S["w_gamma"]
    # @test x_gamma_leja ≈ S["x_gamma_leja"]
    # @test w_gamma_leja ≈ S["w_gamma_leja"]
    # @test x_beta ≈ S["x_beta"]
    # @test w_beta ≈ S["w_beta"]
    # @test x_beta_leja ≈ S["x_beta_leja"]
    # @test w_beta_leja ≈ S["w_beta_leja"]
    # @test x_beta_sym_leja ≈ S["x_beta_sym_leja"]
    # @test w_beta_sym_leja ≈ S["w_beta_sym_leja"]
    # @test x_triang ≈ S["x_triang"]
    # @test w_triang ≈ S["w_triang"]
    end

    @testset "Multi-Index Set Generation" begin
        jj = [2, 3]; min_idx = 0;
        multi_idx_box = create_box_miset(length(jj), jj.+1);

        N = 2;  

        w = 3; 
        multi_idx_TP       = create_tensor_miset(N,w) 
        multi_idx_TD       = create_totaldegree_miset(N,w); 
        # multi_idx_HC       = multiidx_gen(N,rule_HC,w,base); 
        multi_idx_SM       = create_smolyak_miset(N,w);
        rates = [2,3];
        rule_anistropic(l) = sum(rates.*l)
        multi_idx_SM_aniso = create_rule_miset(N,w,rule_anistropic);

        # % fast TD
        #multi_idx_TD_fast = fast_TD_set(N,w); 

        S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_multiidx_set.mat")

        @test get_mi(multi_idx_box) ≈ [S["multi_idx_box"][i,:] .+ fill(1,N) for i in 1:size(S["multi_idx_box"],1)]
        @test get_mi(multi_idx_TP) ≈ [S["multi_idx_TP"][i,:] for i in 1:size(S["multi_idx_TP"],1)]
        @test get_mi(multi_idx_TD) ≈ [S["multi_idx_TD"][i,:] for i in 1:size(S["multi_idx_TD"],1)]
        # @test multi_idx_HC ≈ S["multi_idx_HC"]
        @test get_mi(multi_idx_SM) ≈ [S["multi_idx_SM"][i,:] for i in 1:size(S["multi_idx_SM"],1)]
        # @test multi_idx_SM_aniso ≈ S["multi_idx_SM_aniso"]
        # @test multi_idx_TD_fast ≈ S["multi_idx_TD_fast"]
    end 


# %% test polynomials (tools/polynomials_functions)

# clear

# k = 5; kk = [3,5]; % polynomial order (univariate and bi-variate case, respectively)
# load('test_unit_polynomials','x_unif','X_unif_multidim','x_norm','X_norm_multidim',...
#      'x_exp','X_exp_multidim','x_gamma','X_gamma_multidim','x_beta','X_beta_multidim',...
#      'x_interp','x_eval'); 

# % Legendre and Chebyshev
# a = -2; b = 1; 
# lege_vals = lege_eval(x_unif,k,a,b);
# cheb_vals = cheb_eval(x_unif,k,a,b);

# aa = [-2,1]; bb = [1,2]; 
# multidim_lege_vals = lege_eval_multidim(X_unif_multidim,kk,aa,bb); 
# multidim_cheb_vals = cheb_eval_multidim(X_unif_multidim,kk,aa,bb);

# % Hermite
# mi = 0; sigma = 1;  
# herm_vals = herm_eval(x_norm,k,mi,sigma); 

# mmi = [0,1]; ssigma = [1,0.5]; 
# multidim_herm_vals = herm_eval_multidim(X_norm_multidim,kk,mmi,ssigma); 

# % Laguerre
# lambda = 2; 
# lagu_vals = lagu_eval(x_exp,k,lambda); 

# llambda = [1,2]; 
# multidim_lagu_vals = lagu_eval_multidim(X_exp_multidim,kk,llambda);

# % Generalized Laguerre
# alpha = 1; beta = 2; 
# generalized_lagu_vals = generalized_lagu_eval(x_gamma,k,alpha,beta); % parameters alpha-1,1/beta for consistency with Matlab definitions

# aalpha = [1,2]; bbeta = [1,0.5];  
# multidim_generalized_lagu_vals = generalized_lagu_eval_multidim(X_gamma_multidim,kk,aalpha,bbeta);

# % Jacobi
# a = -2; b = 1; 
# alpha = -0.5; beta = 0.5; 
# jac_vals = jacobi_prob_eval(x_beta,k,alpha,beta,a,b); 

# aa = [-2,1]; bb = [1,2]; 
# aalpha = [-0.5,1]; bbeta = [-0.5,2];   
# multidim_jac_vals = jacobi_prob_eval_multidim(X_beta_multidim,kk,aalpha,bbeta,aa,bb); 

# % univariate interpolant - Lagrange basis
# f = @(x) sin(x); 
# f_interp = f(x_interp); 
# lagr_vals = univariate_interpolant(x_interp,f_interp,x_eval);  

# disp('== testing polynomials ==')
# L = struct('lege_vals',lege_vals,'cheb_vals',cheb_vals,'multidim_lege_vals',multidim_lege_vals,'multidim_cheb_vals',multidim_cheb_vals,...
#            'herm_vals',herm_vals,'multidim_herm_vals',multidim_herm_vals,...
#            'lagu_vals',lagu_vals,'multidim_lagu_vals',multidim_lagu_vals,...
#            'generalized_lagu_vals',generalized_lagu_vals,'multidim_generalized_lagu_vals',multidim_generalized_lagu_vals,...
#            'jac_vals',jac_vals,'multidim_jac_vals',multidim_jac_vals,...
#            'lagr_vals',lagr_vals);
# S = load('test_unit_polynomials',...
#          'lege_vals','cheb_vals','multidim_lege_vals','multidim_cheb_vals',...
#          'herm_vals','multidim_herm_vals',...
#          'lagu_vals','multidim_lagu_vals',...
#          'generalized_lagu_vals','multidim_generalized_lagu_vals',...
#          'jac_vals','multidim_jac_vals', ...
#          'lagr_vals');   
# if isequal_sgmk(L,S)
#     disp('test on polynomials passed')
# end    

# %% test sparse grid generation and reduction (main)

# clear

    @testset "Sparse Grid Generation" begin
        I = [
            [1, 1],
            [1, 2],
            [2, 1],
            [3, 1]
        ];
        I = MISet(I)
        knots1           = GaussLegendrePoints([0,1]);
        knots2           = LejaPoints([-1,1],false,:classic, z->sqrt(0.5));
        lev2knots        = Linear();
        S_given_multiidx = create_sparsegrid(I,knots=[knots1,knots2],rule=[lev2knots,lev2knots]);


        N=2; w=3;
        knots = CCPoints([-1,1]);
        I_smolyak = create_smolyak_miset(N,w);
        S_smolyak = create_sparsegrid(I_smolyak,knots=knots,rule=Doubling());

        # adding one multi-index to S_smolyak
        new_idx = [5, 1];
        I_add = add_mi(I_smolyak,new_idx);
        S_add = create_sparsegrid(I_add,knots=knots,rule=Doubling());
        
        
        # quick preset
        N = 2; w = 3;
        I_quick = create_smolyak_miset(N,w)
        S_quick = create_sparsegrid(I_quick)
        
        S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_grid_gen_and_red.mat")

        compare_sg_to_sgmk(S_given_multiidx, S["S_given_multiidx"])
        compare_sg_to_sgmk(S_smolyak, S["S_smolyak"], S["Sr_smolyak"])
        compare_sg_to_sgmk(S_add, S["S_add"])
        compare_sg_to_sgmk(S_quick, S["S_quick"], S["Sr_quick"])
    end     


    @testset "Sparse Grid Evaluation" begin
        f(x) = sum(x);
        N = 2; w = 3;
        miset = create_smolyak_miset(N,w)
        knots = GaussLegendrePoints([-1,1])
        rule = Linear()
        S  = create_sparsegrid(miset,knots=knots,rule=rule);

        f_evals = f.(sort!(get_grid_points(S)));

        w = 4;
        miset = create_smolyak_miset(N,w)
        T = create_sparsegrid(miset,knots=knots,rule=rule);
        f_evals_rec = f.(sort!(get_grid_points(T)));

        S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_evaluate.mat")
        @test f_evals ≈ vec(S["f_evals"])
        @test f_evals_rec ≈ vec(S["f_evals_rec"])
    end

    @testset "Sparse Grid Integration" begin
        f(x) = prod(1.0 ./sqrt.(x.+3));

        N = 4; w =4 ;
        knots = CCPoints([-1,1]);
        rule = Doubling()
        miset = create_smolyak_miset(N,w)

        S = create_sparsegrid(miset,knots=knots,rule=rule);

        f_quad = integrate_on_sparsegrid(S,f.(get_grid_points(S)));

        S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_quadrature.mat")
        @test f_quad ≈ S["f_quad"]
    end

    @testset "Sparse Grid Interpolation" begin
        f(x) = prod(1.0 ./sqrt.(x.+3));

        N = 2; w = 4;
        knots = CCPoints([-1,1]);
        rule = Doubling()
        miset = create_smolyak_miset(N,w)

        S = create_sparsegrid(miset,knots=knots,rule=rule);

        x = range(-1, stop=1, length=10)
        non_grid_points = [[x1,x2] for x1 in x, x2 in x]
        f_on_grid = f.(get_grid_points(S));
        f_values = interpolate_on_sparsegrid(S,f_on_grid,non_grid_points);

        S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_interpolate.mat")
        @test vec(f_values) ≈ vec(S["f_values"])
    end

    @testset "Polynomial Chaos Expansions" begin
        f(x) = prod(1.0 ./sqrt.(x.+3));

        N = 2; w = 5; a = -1; b = 1;
        knots = CCPoints([a,b]);
        rule = Linear()
        idxset(mi) = prod(mi)
        miset = create_rule_miset(N,w,idxset)
        S = create_sparsegrid(miset,knots=knots,rule=rule);

        f_on_grid = f.(get_grid_points(S));

        sga = SparseGridApproximation(S, f_on_grid);
        sgs = convert_to_spectral_approximation(sga);

        S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_gPCE.mat")
        
        # @debug "Skipping PCE comparison, we only use Chebyshev polynomials at the moment not Legendre"
        # @test sg ≈ S["PCE_coeffs"]
    end

    # %% test on Sobol indices 

    # clear

    # f = @(x) 1./(1 + 5*x(1,:).^2 + x(2,:).^2 + x(3,:).^2); 

    # N = 3; w = 5; a = -1; b = 1; 
    # knots     = @(n) knots_uniform(n,a,b); 
    # lev2knots = @lev2knots_lin; 
    # idxset    = @(i) prod(i); 
    # S         = create_sparse_grid(N,w,knots,lev2knots,idxset);
    # Sr        = reduce_sparse_grid(S);
    # f_on_grid = evaluate_on_sparse_grid(f,Sr);

    # [Sob_i,Tot_Sob_i,Mean,Var] = compute_sobol_indices_from_sparse_grid(S,Sr,f_on_grid,[a a a; b b b],'legendre');

    # disp('== testing Sobol indices ==')
    # L = struct('Sob_i',Sob_i,'Tot_Sob_i',Tot_Sob_i,'Mean',Mean,'Var',Var);
    # S = load('test_unit_sobol','Sob_i','Tot_Sob_i','Mean','Var');   
    # if isequal_sgmk(L,S)
    #     disp('test on Sobol indices passed')
    # end    

    # %% test on gradient and Hessian

    # clear

    # f = @(x) 1./(1+0.5*sum(x.^2)); 

    # N = 2; aa = [4 1]; bb = [6 5]; domain = [aa; bb]; 
    # w = 4;
    # knots1    = @(n) knots_CC(n,aa(1),bb(1));
    # knots2    = @(n) knots_CC(n,aa(2),bb(2));
    # S         = create_sparse_grid(N,w,{knots1,knots2},@lev2knots_doubling);
    # Sr        = reduce_sparse_grid(S);
    # f_on_grid = evaluate_on_sparse_grid(f,Sr);

    # x1        = linspace(aa(1),bb(1),10);
    # x2        = linspace(aa(2),bb(2),30);
    # [X1,X2]   = meshgrid(x1,x2); 
    # eval_points = [X1(:),X2(:)]'; 

    # grad = derive_sparse_grid(S,Sr,f_on_grid,domain,eval_points);
    # eval_point_hess = eval_points(:,10); % pick one point
    # Hess = hessian_sparse_grid(S,Sr,f_on_grid,domain,eval_point_hess);

    # disp('== testing gradients and derivatives ==')
    # L = struct('grad',grad,'Hess',Hess);
    # S = load('test_unit_gradient_and_hessian','grad','Hess');   
    # if isequal_sgmk(L,S)
    #     disp('test on gradients and derivatives passed')
    # end    

    @testset "Adaptive Sparse Grid Approximation" begin

        f(x) = 1.0 ./ (x[1].^2 .+x[2].^2 .+ 0.3);

        N = 2; a = -1; b = 1;
        knots     = CCPoints([a,b]);
        lev2knots = Doubling();
        controls = Dict(
            :maxpts => 200,
            :proftol => 1e-10,
        )

        (S_Linf, f_on_z) = adaptive_sparsegrid(f,N; maxpts=controls[:maxpts], proftol=controls[:proftol], rule=lev2knots,knots=knots, type=:Linf);
        (S_Linfcost, f_on_z) = adaptive_sparsegrid(f,N; maxpts=controls[:maxpts], proftol=controls[:proftol] ,rule=lev2knots,knots=knots, type=:Linfcost);
        (S_deltaint, f_on_z) = adaptive_sparsegrid(f,N; maxpts=controls[:maxpts], proftol=controls[:proftol], rule=lev2knots,knots=knots, type=:deltaint);
        (S_deltaintcost, f_on_z) = adaptive_sparsegrid(f,N; maxpts=controls[:maxpts], proftol=controls[:proftol], rule=lev2knots,knots=knots, type=:deltaintcost);

        # N = 3; 
        # f = @(x) prod(x);
        # knots     = @(n) knots_CC(n, 0, 2);
        # lev2knots = @lev2knots_doubling;
        # controls.nested = true;

        # controls.paral     = NaN; 
        # controls.max_pts   = 200;
        # controls.prof_toll = 1e-10;
        # prev_adapt         = [];
        # controls.plot      = false;

        # controls.var_buffer_size = 2;
        # S_buff = adapt_sparse_grid(f,N,knots, lev2knots, prev_adapt, controls);

        S = read_sgmk_mat("https://raw.githubusercontent.com/lorenzo-tamellini/sparse-grids-matlab-kit/main/docs-examples/testing_unit/test_unit_adaptive.mat")

        # compare_sg_to_sgmk(S_Linf, S["S_Linf"]["S"], S["S_Linf"]["Sr"])
        # compare_sg_to_sgmk(S_Linfcost, S["S_Linf_newpts"]["S"], S["S_Linf_newpts"]["Sr"])
        # compare_sg_to_sgmk(S_deltaint, S["S_deltaint"]["S"], S["S_deltaint"]["Sr"])
        # compare_sg_to_sgmk(S_deltaintcost, S["S_deltaint_newpts"]["S"], S["S_deltaint_newpts"]["Sr"])

        # We cannot directly compared to SGMK as the algorithms are different - the reduced margin contitutes part of the approximation in the whilst in SparseGridsKit it does not.
        # We verify the adaptive approximation is accurate wrt the true function
        test_points = [[x1,x2] for x1 in range(a, stop=b, length=10^2) for x2 in range(a, stop=b, length=10^2)]
        f_test = f.(test_points)
        f_Linf = interpolate_on_sparsegrid(S_Linf, f_on_z, test_points)
        f_Linfcost = interpolate_on_sparsegrid(S_Linfcost, f_on_z, test_points)
        f_deltaint = interpolate_on_sparsegrid(S_deltaint, f_on_z, test_points)
        f_deltaintcost = interpolate_on_sparsegrid(S_deltaintcost, f_on_z, test_points)

        tol = 1e-2;
        @test all(isapprox(f_test[i], f_Linfcost[i]; atol=tol) for i in 1:length(f_test))
        @test all(isapprox(f_test[i], f_Linf[i]; atol=tol) for i in 1:length(f_test))
        @test all(isapprox(f_test[i], f_deltaint[i]; atol=tol) for i in 1:length(f_test))
        @test all(isapprox(f_test[i], f_deltaintcost[i]; atol=tol) for i in 1:length(f_test))
    end
end