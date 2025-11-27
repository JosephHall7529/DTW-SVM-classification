using MLJModelInterface
using MLJBase

abstract type DTW_results end

struct DTW_train_result <: DTW_results
    metadata::DataFrame
    cosine_cost::DataFrame
    flat_top_cost::DataFrame
end

struct DTW_test_result <: DTW_results
    metadata::DataFrame
    cosine_cost::DataFrame
    flat_top_cost::DataFrame
end

mutable struct DTW_data
    tok_shots::Vector{Tuple{String, Int}}
    features::Vector{String}
    short_IP::Bool
    shot_ind_dict::Dict{Tuple, Int}
    ind_shot_dict::Dict{Int, Tuple}
    target_shot_length::Int
    profile_data::Dict{Tuple, Array}
    flat_top_data::Dict{Tuple, Array}
    training_data::DTW_results
    testing_data::DTW_results
    function DTW_data(features::Vector{String})
        # require IP and PNBI for the two comparable durations 
        features_ = vcat(features, ["IP", "PNBI"])
        features_ = unique(features_)

        tss = tok_shots((@subset which(features_) @byrow in(:tok, ["aug"])))

        ℓT = length(tss)
        ind_shot_dict = Dict(1:ℓT .=> tss)
        shot_ind_dict = Dict(tss .=> 1:ℓT)

        new(tss, sort!(features), false, shot_ind_dict, ind_shot_dict, 0)
    end
end

mutable struct DTW_SVM <: MLJBase.Deterministic
    target_shot_length::Int
    transportcost::Float64
    C_cosine::Float64
    C_flat_top::Float64
    C_cost::Float64
    short_IP::Bool
    database::DTW_data
end
function DTW_SVM(features::Vector{String};
                    target_shot_length=100,
                    transportcost=1.0,
                    C_cosine=1.0,
                    C_flat_top=1.0,
                    C_cost=1.0,
                    short_IP=false)
    database = DTW_data(features)
    DTW_SVM(target_shot_length, transportcost, C_cosine, C_flat_top, C_cost, short_IP, database)
end

function least_data_length(data_::DTW_data)
    dict = Dict{Tuple{String, Int}, Dict{String, NamedTuple}}() 
    if data_.short_IP
        for ts in data_.tok_shots
            IP_t_size = [sum(FTIPs[ts][1] .< profiles(ts..., feat).t .< FTIPs[ts][2]) for feat in data_.features]
            IP_smallest_size = minimum(IP_t_size)
            NBI_t_size = [sum(FTNBI[ts][1] .< profiles(ts..., feat).t .< FTNBI[ts][2]) for feat in data_.features]
            NBI_smallest_size = minimum(NBI_t_size)
            dict[ts] = Dict("IP" => (feat = data_.features[IP_t_size .== IP_smallest_size][1], length = IP_smallest_size),
                            "NBI" => (feat = X.features[NBI_t_size .== NBI_smallest_size][1], length = NBI_smallest_size)
            )
        end
    else
        for ts in data_.tok_shots
            IP_t_size = [sum(FTIP[ts][1] .< profiles(ts..., feat).t .< FTIP[ts][2]) for feat in data_.features]
            IP_smallest_size = minimum(IP_t_size)
            NBI_t_size = [sum(FTNBI[ts][1] .< profiles(ts..., feat).t .< FTNBI[ts][2]) for feat in data_.features]
            NBI_smallest_size = minimum(NBI_t_size)
            dict[ts] = Dict("IP" => (feat = data_.features[IP_t_size .== IP_smallest_size][1], length = IP_smallest_size),
                            "NBI" => (feat = data_.features[NBI_t_size .== NBI_smallest_size][1], length = NBI_smallest_size)
            )
        end
    end
    return dict
end

function fill_dtw_data!(data_::DTW_data, target_length::Int)
    ℓf = length(data_.features)
    LSD = least_data_length(data_)

    profile_data = Dict{Tuple, Array{Float64}}()
    flat_top_data = Dict{Tuple, Array{Float64}}()

    for (i, ts) in enumerate(data_.tok_shots)
        Random.seed!(i)
        mat_ip = Float64[]
        mat_nbi = Float64[]
        IP_step_constraint = LSD[ts]["IP"] 
        if IP_step_constraint.length > target_length
            IP_step_constraint = merge(IP_step_constraint, (length=target_length,))
        end
        NBI_step_constraint = LSD[ts]["NBI"] 
        if NBI_step_constraint.length > target_length
            NBI_step_constraint = merge(NBI_step_constraint, (length=target_length,))
        end

        constraint_feature_IP = profiles(ts..., IP_step_constraint.feat).t
        if data_.short_IP
            constraint_time_IP = constraint_feature_IP[FTIPs[ts][1] .< constraint_feature_IP .< FTIPs[ts][2]]
        else
            constraint_time_IP = constraint_feature_IP[FTIP[ts][1] .< constraint_feature_IP .< FTIP[ts][2]]
        end
        constraint_step_IP = div(length(constraint_time_IP), IP_step_constraint.length)
        constraint_time_IP = constraint_time_IP[1:constraint_step_IP:end] 
        constraint_feature_NBI = profiles(ts..., NBI_step_constraint[1]).t
        constraint_time_NBI = constraint_feature_NBI[FTNBI[ts][1] .< constraint_feature_NBI .< FTNBI[ts][2]]
        constraint_step_NBI = div(length(constraint_time_NBI), NBI_step_constraint.length)
        constraint_time_NBI = constraint_time_NBI[1:constraint_step_NBI:end] 
        for (n, feat) in enumerate(data_.features)
            P = profiles(ts..., feat)

            ind = Vector{Int}()
            for t in constraint_time_IP
                ind_last = findlast(i -> i <= t, P.t)
                if ind_last == nothing
                    push!(ind, findfirst(i -> i >= t, P.t))
                else
                    push!(ind, ind_last)
                end
            end
            append!(mat_ip, abs.(P.y[ind, 1].*normalise_2D_features[feat]))

            ind = Vector{Int}()
            for t in constraint_time_NBI
                ind_last = findlast(i -> i <= t, P.t)
                if ind_last == nothing
                    push!(ind, findfirst(i -> i >= t, P.t))
                else
                    push!(ind, ind_last)
                end
            end
            append!(mat_nbi, abs.(P.y[ind, 1].*normalise_2D_features[feat]))

            if n == ℓf
                append!(mat_ip, constraint_time_IP) 
                append!(mat_nbi, constraint_time_NBI)
            end
        end
            
        profile_data[ts] = reshape(mat_ip, :, ℓf+1)'
        flat_top_data[ts] = reshape(mat_nbi, :, ℓf+1)'
    end
    data_.profile_data = profile_data
    data_.flat_top_data = flat_top_data

    return data_
end

function model_fill_train!(model::DTW_SVM, X, y)
    ℓT = length(X) #cost matrices are equal height and width of size ℓT
    indices = [model.database.shot_ind_dict[ts] for ts in X]
    @assert length(indices) == ℓT "One of the training shots cannot be located in the dictionary shot -> ind"

    cosine_cost = zeros(ℓT, ℓT)
    flat_top_cost = zeros(ℓT, ℓT)
    
    # find indices of training data, initialise a dataframe with the metadata for training
    metadata = DataFrame(:order => 1:ℓT, :ind => indices, :tok_shots => X, :y => y)
    sort!(metadata, [:y, :tok_shots], rev=[true, false])

    path_cos_ℓ = zeros(ℓT, ℓT)
    path_ft_ℓ = zeros(ℓT, ℓT)
    for (i, train_ts_i) in ProgressBar(enumerate(metadata.tok_shots))
        for (j, train_ts_j) in ProgressBar(enumerate(metadata.tok_shots))
            data_j = model.database.profile_data[train_ts_j]
            data_i = model.database.profile_data[train_ts_i]

            COST_cos, path_cos_a, path_cos_b  = dtw(data_j, data_i, CosineDist(); transportcost=model.transportcost)
            cosine_cost[i, j] = COST_cos
            path_cos_ℓ[i, j] = length(path_cos_a)

            data_j = model.database.flat_top_data[train_ts_j]
            data_i = model.database.flat_top_data[train_ts_i]
        
            COST_ft, path_ft_a, path_ft_b = dtw(data_j, data_i, Euclidean(); transportcost=model.transportcost)
            flat_top_cost[i, j] = COST_ft
            path_ft_ℓ[i, j] = length(path_ft_a)
        end
    end

    cosine_cost = round.(cosine_cost ./ path_cos_ℓ, sigdigits=5)
    flat_top_cost = round.(flat_top_cost ./ path_ft_ℓ, sigdigits=5) 
    
    df_names = df_ts_naming(metadata.tok_shots)

    cosine_cost = DataFrame(cosine_cost[:, :], df_names)
    flat_top_cost = DataFrame(flat_top_cost[:, :], df_names) 

    insertcols!(cosine_cost, 1, :shots => metadata.tok_shots) 
    insertcols!(flat_top_cost, 1, :shots => metadata.tok_shots)

    model.database.training_data = DTW_train_result(metadata, cosine_cost, flat_top_cost)
    return model.database
end
function model_fill_test!(model::DTW_SVM, X)
    training = model.database.training_data
    ℓtr = size(training.metadata, 1)
    ℓT = length(X) #cost matrices are equal height and width of size ℓT
    indices = [model.database.shot_ind_dict[ts] for ts in X]
    @assert length(indices) == ℓT "One of the training shots cannot be located in the dictionary shot -> ind"

    cosine_cost = zeros(ℓtr, ℓT)
    flat_top_cost = zeros(ℓtr, ℓT)
    
    # find indices of training data, initialise a dataframe with the metadata for training 
    metadata = DataFrame(:order => 1:ℓT, :ind => indices, :tok_shots => X)

    path_cos_ℓ = zeros(ℓtr, ℓT) 
    path_ft_ℓ = zeros(ℓtr, ℓT) 
    for (i, train_ts) in ProgressBar(enumerate(training.metadata.tok_shots))
        for (j, test_ts) in enumerate(metadata.tok_shots)
            data_j = model.database.profile_data[test_ts]
            data_i = model.database.profile_data[train_ts]

            COST_cos, path_cos_a, path_cos_b  = dtw(data_j, data_i, CosineDist(); transportcost=model.transportcost)
            cosine_cost[i, j] = COST_cos
            path_cos_ℓ[i, j] = length(path_cos_a)
            
            data_j = model.database.flat_top_data[test_ts]
            data_i = model.database.flat_top_data[train_ts]

            COST_ft, path_ft_a, path_ft_b = dtw(data_j, data_i, Euclidean(); transportcost=model.transportcost)
            flat_top_cost[i, j] = COST_ft
            path_ft_ℓ[i, j] = length(path_ft_a)
        end
    end

    cosine_cost = round.(cosine_cost ./ path_cos_ℓ, sigdigits=5)
    flat_top_cost = round.(flat_top_cost ./ path_ft_ℓ, sigdigits=5) 
    
    df_names = df_ts_naming(metadata.tok_shots)

    cosine_cost = DataFrame(cosine_cost[:, :], df_names)
    flat_top_cost = DataFrame(flat_top_cost[:, :], df_names) 

    insertcols!(cosine_cost, 1, :train => training.metadata.tok_shots) 
    insertcols!(flat_top_cost, 1, :train => training.metadata.tok_shots)

    model.database.testing_data = DTW_test_result(metadata, cosine_cost, flat_top_cost)
    return model
end

dtw_svm_report(fitresult) = (SVs = fitresult.SVs.nSV)

function MLJBase.fit(model::DTW_SVM, verbosity::Int, X, y) 
    if model.database.target_shot_length !== model.target_shot_length
        fill_dtw_data!(model.database, model.target_shot_length)
    end
    model_fill_train!(model, X, y)
    Xmatrix = (ftX = Array(model.database.training_data.flat_top_cost[!, 2:end]),
                cosX = Array(model.database.training_data.cosine_cost[!, 2:end]))

    K = exp.(-(Xmatrix.cosX ./ (model.C_cosine)).^2) .*
        exp.(-(Xmatrix.ftX ./ (model.C_flat_top)).^2)

    X = MLJBase.matrix(K)
    y = [y[i] for i in model.database.training_data.metadata.order]
    fitresult = svmtrain(X, y, kernel=Kernel.Precomputed, cost=model.C_cost)

    report = dtw_svm_report(fitresult)
    # report = nothing
    cache = nothing

    return fitresult, cache, report
end

function MLJBase.predict(model_::DTW_SVM, fitresult, Xnew)
    if model_.database.target_shot_length !== model_.target_shot_length
        fill_dtw_data!(model_.database, model_.target_shot_length)
    end
    model_fill_test!(model_, Xnew)
    test_data = model_.database.testing_data.metadata

    Xmatrix = (ftX = Array(model_.database.testing_data.flat_top_cost[!, 2:end]), 
                cosX = Array(model_.database.testing_data.cosine_cost[!, 2:end]))

    println(size(Xmatrix.ftX))
    K = exp.(-(Xmatrix.cosX ./ (model_.C_cosine)).^2) .* 
        exp.(-(Xmatrix.ftX ./ (model_.C_flat_top)).^2) 

    X = MLJBase.matrix(K)
    ỹ, confidence = svmpredict(fitresult, X)
    test_data.y = [ỹ[i] for i in test_data.order]
    test_data.confidence = [confidence[:, i] for i in test_data.order]

    sort!(test_data, [:y, :tok_shots], rev=[true, false])

    return ỹ
end

function half_model(features::Vector{String}, X, y;
                target_shot_length=100,
                transportcost=1.0,
                C_cosine=1.0,
                C_flat_top=1.0,
                C_cost=1.0,
                short_IP=false,
                kernel::Bool=false)
    model = DTW_SVM(features;
                target_shot_length=target_shot_length,
                transportcost=transportcost,
                C_cosine=C_cosine,
                C_flat_top=C_flat_top,
                C_cost=C_cost,
                short_IP=short_IP)
    fill_dtw_data!(model.database, target_shot_length)
    model_fill_train!(model, X, y)

    if kernel
        Xmatrix = (ftX = Array(model.database.training_data.flat_top_cost[! , 2:end]),
                cosX = Array(model.database.training_data.cosine_cost[!, 2:end]))

        K = exp.(-(Xmatrix.cosX ./ (model.C_cosine)).^2) .*
            exp.(-(Xmatrix.ftX ./ (model.C_flat_top)).^2)

        return K, model
    end
    return model
end
function C_extremes(features::Vector{String}, target_length::Int, X, y)
    
    ℓ = length(X)

    model = half_model(features, target_length, X, y)
    CC = Array(model.database.training_data.cosine_cost[:, 2:end])
    FT = Array(model.database.training_data.flat_top_cost[:, 2:end])

    results = Dict()
    for (df, name) in zip([CC, FT], ["cosine_cost", "flat_top_cost"])
        results[name] = extrema(vcat([[df[i, j] for j in i+1:ℓ] for i in 1:ℓ-1]...))
    end
    return results
end
function C_extremes(model::DTW_SVM)
    CC = Array(model.database.training_data.cosine_cost[:, 2:end])
    FT = Array(model.database.training_data.flat_top_cost[:, 2:end])

    ℓ = size(CC, 1)

    results = Dict()
    for (df, name) in zip([CC, FT], ["cosine_cost", "flat_top_cost"])
        results[name] = extrema(vcat([[df[i, j] for j in i+1:ℓ] for i in 1:ℓ-1]...))
    end
    return results
end