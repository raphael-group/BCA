include("Utils.jl")
using DataFrames, Serialization, ManifoldLearning

function get_sim(n, d, sampling)
    A_cycle, A_type = 10, 10 # 10, 10 and 10, 4
    X, cycle_labs, type_labs, t = simple_sim(n, A_cycle, A_type, σ_cycle=0.1, σ_type=0.3, use_seed=false, sampling=sampling); 
    σ_noise = 5
    X = hcat(X, rand(Normal(0, σ_noise), n, d));
    return X, cycle_labs, type_labs, t 
end

function main(nsims, n, d, sample_type)
    Random.seed!(123)
    ss_true, ss_pc, ss_ours, ss_lda, ss_lem = [], [], [], [], []
    for i in 1:nsims
        num_pcs = 3

        X, cycle_labs, type_labs, t = get_sim(n, d, :normal)
        pcs = collect(getPCs(collect(X'), true)');
        Y, R = cluster_projection(collect(X'),cycle_labs, false, :mean);
        mlda, Y_lda = mLDA(collect(X'), cycle_labs, true);
        M = fit(LEM, collect(pcs[:,1:num_pcs]'), maxoutdim=2)
        R = ManifoldLearning.transform(M)

        s_true = slingshot_ordering(collect(X[:,1:2]'), cycle_labs)
        s_pc = slingshot_ordering(collect(pcs[:,1:num_pcs]'), cycle_labs)
        s_ours = slingshot_ordering(Y, cycle_labs)
        s_lda = slingshot_ordering(Y_lda, cycle_labs);
        s_lem = slingshot_ordering(R, cycle_labs);

        push!(ss_true, fix_angle(s_true, t)[1])
        push!(ss_pc, fix_angle(s_pc, t)[1])
        push!(ss_ours, fix_angle(s_ours, t)[1])
        push!(ss_lda, fix_angle(s_lda, t)[1])
        push!(ss_lem, fix_angle(s_lem, t)[1])
    end
    labels = repeat(["True", "PCA", "Ours", "mLDA", "LEM"], inner=nsims)
    results = vcat(ss_true, ss_pc, ss_ours, ss_lda, ss_lem)
    df = DataFrame(Methods = labels, KendallsTau=results)
    serialize("taus-$(nsims)-$(n)-$(d)-uniform.jls", df)
end

main(100, 1_200, 1_000, :uniform)