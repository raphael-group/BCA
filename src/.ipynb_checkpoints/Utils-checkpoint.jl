#import Pkg; Pkg.activate("../."); #Pkg.instantiate()
using LinearAlgebra, MultivariateStats, Statistics, DataFrames, CSV, StatsBase, Random, RCall

# A useful trick to transpose a dataframe
transposedf(df) = DataFrame([[names(df)]; collect.(eachrow(df))], [:column; Symbol.(axes(df, 1))])

# Makes the columns of the matrix X have empircal mean of 0
meancenter(X, dim=2) = X .- mean(X, dims=dim)

# Scale columns of X to have mean 0 and var 1, i.e. Z-score them. 
function scale_data(X, dim=1)
    dt = StatsBase.fit(StatsBase.ZScoreTransform, X, dims=dim)
    return  StatsBase.transform(dt, X)
end

# Takes a matrix X that is dim reduced that is reduced dim × samples, transposes it and makes it into a Dataframe that is of size samples × reduced dim
to_dataframe(X, colnames=["PC_$(i)" for i in 1:size(X,1)]) = DataFrame(collect(X'), colnames)

# TODO: add type colors
function addlabels!(df, cycle_labs, type_labs, cycle_colors=(:black, :yellow, :red))
    @assert cycle_labs isa Vector{Int} "Cycle labels is not a vector of ints."
    @assert type_labs isa Vector{Int} "Type labels is not a vector of ints."
    
    K_cycle = length(unique(cycle_labs))
    K_type = length(unique(type_labs))
    cycle_to_colors = Dict(i-1 => cycle_colors[i] for i in 1:K_cycle)
    clabs = [cycle_to_colors[l] for l in cycle_labs]
    tlabs = type_labs # TODO: add dictionary for more fine grained dictionary
    insertcols!(df, 1, :cell_cycle => clabs)
    insertcols!(df, 1, :cell_type => tlabs)
end

# Input: `X` genes by cells matrix, `labels` vector of cell cluster labels
# Output: `Y` cells by |`labels`| matrix either as a dataframe or a matrix depending on value of `as_df`
function MBCVE(X, labels; weighted=false, as_df=false)
    X = meancenter(X)
    df = to_dataframe(X)
    insertcols!(df, 1, :labels => labels)
    gdf = groupby(df, :labels)
    if weighted
        M = combine(gdf, valuecols(gdf) .=> mean)
        clust_sizes = diagm([size(g,1) for g in gdf]) # |C_1|, \dots, |C_K|
        M_tilde = clust_sizes * Matrix(M[:,2:end]) # without labels = 2:end
        U_m, S_m, V_m = svd(M_tilde)
        Y = V_m' * X
    else
        M = combine(gdf, valuecols(gdf) .=> mean)
        U_m, S_m, V_m = svd(Matrix(M[:,2:end]))
        Y = V_m' * X
    end
    if as_df
        return insertcols!(to_dataframe(Y, ["BCC_$(i)" for i in 1:size(Y,1)]), 1, :labels => labels)
    else
        return Y, V_m
    end
end


# Assume X is genes x cells
# Returns PCs as a dataframe
function getPCs(X, return_as_mat=false)
    PCA_model = MultivariateStats.fit(PCA, X, pratio=0.9) # keep all pcs with 90% of variance
    pcs = MultivariateStats.transform(PCA_model, X)
    if return_as_mat
        return Matrix(pcs)
    end
    df = to_dataframe(pcs)
    return df 
end

# X genes × cells; y integer labels of size cells
# Output: `Y` mLDA embedding of X size cells by lda features.  
function mLDA(X, y, return_S=false)
    model = MulticlassLDA  
    k = length(unique(y))
    label_dict = Dict(l => i for (i,l) in enumerate(unique(y)))
    labels = [label_dict[l] for l in y]
    mlda = MultivariateStats.fit(model, k+1, X, labels)
    Y = MultivariateStats.transform(mlda, X)
    if !return_S
        return DataFrame(collect(Y'), :auto)
    else
       return mlda, Matrix(Y)
    end
end

###################### Simulation Utilities #########################
# Compute column variances of X. 
get_variances(X) = [var(X[:,i]) for i in 1:size(X,2)]

# Computese unweighted between cluster variance of clusters found in labs with respect to X. 
# If X is meancentered, then this is equivalent to eq (10) in the paper, since x̄ = zero vector. 
function bclust_var(X, labs)
    A = cluster_indicator(labs, false).^2
    if size(X,2) > 1
        M = A' * meancenter(X) # k × features
    else
        M = A' * (X .- mean(X))
    end
    k = size(M,1)
    return 1/(k-1) * sum(norm(M[i,:])^2 for i in 1:k)
end

# Computes within cluster variance of clusters found in labs with respect to X. 
# Assumes X is of samples by features
function wclust_var(X, labs)
      return total_var(X) - bclust_var(X, labs)                              
end 

# Computes total variance of X. 
# Assumes X is of samples by features
total_var(X) = sum(var(x) for x in eachcol(X))  

# A utility function that partitions a vector t \in [0,hi] into num_clusters equally
function get_cycle_labels(t, num_clusters, hi)
    n = length(t) 
    labels = zeros(Int, n)  
    for i in 1:num_clusters 
        a, b = i*hi/num_clusters, (i+1)*hi/num_clusters
        inds = findall(x -> a<= x <= b, t) 
        labels[inds] .= i
    end
    return labels
end

# Our main simulation function. 
function simple_sim(n, A_circle, A_blob; 
        K_cycle=3, σ_cycle=0.1, σ_type=0.3, 
        scale_vars=false, use_seed=true, 
        rng=123, verbose=false, sampling=:uniform)
    if use_seed
        Random.seed!(rng)
    end
    maxval = 2*π
    if sampling == :uniform
        t = rand(Uniform(0,maxval),n) # t[i] ∼ Unif(0,2π)
    elseif sampling == :normal
        sigs = rand(Uniform(0, 2), 3) # sample variances uniformly from [0,2]
        @assert isinteger(n/3) "n needs to be divisible by 3 for normal sampling."
        ns = Int(n/3)
        t = mod.(vcat([rand(Normal(π/4, sigs[1]), ns), rand(Normal(3π/4, sigs[2]), ns), rand(Normal(3π/2, sigs[3]), ns)]...), maxval)
    end
    
    cycle_labels = get_cycle_labels(t, K_cycle, maxval)
    ν_cycle_x = rand(Normal(0,σ_cycle), n) # ν[i] ∼ N(0,var=σ_cycle)
    ν_cycle_y = rand(Normal(0,σ_cycle), n) # ν[i] ∼ N(0,var=σ_cycle)
    x = A_circle * cos.(t) + ν_cycle_x
    y = A_circle * sin.(t) + ν_cycle_y
    
    z_labels = A_blob * (1 .- 2rand(Bernoulli(0.5), n)) # z[i] ∼ Bernoulli(0.5) with support shifted from {0,1} → {-A_blob, A_blob}
    z = z_labels + rand(Normal(0,σ_type), n)
    if scale_vars 
        X = scale_data(meancenter([x y z], 1))
    else
        X = meancenter([x y z], 1)
    end
    if verbose
        println("Var for Cell Type: $(round(bclust_var(X, z_labels), digits=4)); Var for Cell Cycle: $(round(bclust_var(X, cycle_labels), digits=4))")
    end

    return X, cycle_labels, z_labels, t
end

# Iterate over all starting positions and compute maximum kendall's tau
function find_best_score(ts, true_ts, curr_score)
    curr_start = zeros(length(ts))
    for t in ts
        test_ts = mod.(ts .- t, 1.0)
        test_score = corkendall(test_ts, true_ts)
        if test_score > curr_score
            curr_score, curr_start = test_score, t
        end
    end
    return curr_score, curr_start
end

# Compute corrected (for rotation) kendall's tau betweeen inferred ts (pseudotime), and ground truth ts (true_ts). 
function fix_angle(ts, true_ts)
    ts_normed = ts / (maximum(ts) - minimum(ts)) #normlaize to [0,1]
    best_start = ts_normed[1] # initialize best_start
    s1_init = corkendall(ts_normed, true_ts)
    s2_init = corkendall(-1*ts_normed, true_ts)
    score1, start1 = find_best_score(ts_normed, true_ts, s1_init)
    score2, start2 = find_best_score(-1*ts_normed, true_ts, s2_init)
    if (score1 > score2) best_start = start1 else best_start = start2 end 
    return max(score1, score2), mod.(ts_normed .- best_start, 1.0)
end                                                            
 
# Input: X (reduced dims by cells) matrix, cluster_labels: vector of length cells.
# Return: vector of length n of pseudotimes inferred by Slingshot                                                                          
function slingshot_ordering(X, cluster_labels, start_clust=1)
    R"library(slingshot)"
    R"sds <- slingshot($(collect(X')), $(cluster_labels), start.clus=$(start_clust))"
    R"ts <- slingPseudotime(sds)"
    @rget ts
    return vec(ts)
end
                            

#########################################################


# requires StatsBase for countmap; 
# labels: an n dimensional vector that contains cluster labels
# return_labels: boolean that determines whether to return a dictionary about the discrete mapping of cluster labels
# Returns: n×k matrix cluster indicator matrix A s.t. A[i, k] = 1/sqrt(|C_k|) if and only if i ∈ C_k; 0 otherwise. 
function cluster_indicator(labels, return_labs=false)
    d = countmap(labels)    
    k = length(keys(d)) # number of clusters
    n = length(labels)
    A = zeros(n, k)
    cluster_ints = Dict(k => i for (i, k) in enumerate(keys(d))) # remap clusters to ints
    for i in 1:n
        cluster_assignment = labels[i]
        A[i,cluster_ints[cluster_assignment]] = 1 / sqrt(d[cluster_assignment]) 
    end
    if return_labs 
        return A, cluster_ints
    else
        return A
    end
end

