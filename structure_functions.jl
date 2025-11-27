function original_space(data::hybrid_classification)
    return data.data_0D.data.original_space
end
function original_space(data::data0D)
    return data.data.original_space
end
function original_space(data::global_data)
    return data.original_space
end

function log_space(data::hybrid_classification)
    return data.data_0D.data.log_space
end
function log_space(data::data0D)
    return data.data.log_space
end
function log_space(data::global_data)
    return data.log_space
end

function baseline(data::hybrid_classification)
    return data.data_0D.baseline
end
function hybrid(data::hybrid_classification)
    return data.data_0D.hybrid
end
function tbd(data::hybrid_classification)
    return data.data_0D.tbd
end
function baseline(data::data0D)
    return data.baseline
end
function hybrid(data::data0D)
    return data.hybrid
end
function tbd(data::data0D)
    return data.tbd
end

function profiles()
    return DATA_.data_2D.data
end
function profiles(tok::String)
    return DATA_.data_2D.data[tok]
end
function profiles(tok::String, shot::Int64)
    return DATA_.data_2D.data[tok].profiles[shot]
end
function profiles(tok, shot::Int64, feat::String)
    profile_shot = profiles(tok, shot) 
    if typeof(profile_shot) == ITPA_profile_extract
        return profile_shot.data[feat]
    elseif typeof(profile_shot) == SAL_profile_extract
        return profile_shot.data[SAL_dict[feat]][feat]
    elseif typeof(profile_shot) == AUG_profile_extract
        return profile_shot.data[AUG_dict[feat]][feat]
    end 
end
function profiles(data::data2D)
    return data.data
end
function profiles(data::data2D, tok::String)
    return data.data[tok]
end
function profiles(data::data2D, tok::String, shot::Int64) 
    return data.data[tok].profiles[shot]
end
function profiles(data::data2D, tok, shot::Int64, feat::String)
    profile_shot = profiles(data, tok, shot) 
    if typeof(profile_shot) == ITPA_profile_extract
        return profile_shot.data[feat]
    elseif typeof(profile_shot) == SAL_profile_extract
        return profile_shot.data[SAL_dict[feat]][feat]
    elseif typeof(profile_shot) == AUG_profile_extract
        return profile_shot.data[AUG_dict[feat]][feat]
    end 
end


# profile data
function id(ts::Tuple{String, Int})
    id = String((@subset id_codes_2D @byrow begin 
        :TOK == ts[1]
        in(:SHOT, ts[2])
    end).id[1])
end
function id!(data::data2D)
    for row in eachrow(id_codes_2D)
        tok = row.TOK
        shot = row.SHOT
        prof = profiles(data, tok, shot)
        prof.id = String((@subset id_codes_2D @byrow begin 
            :TOK == tok
            in(:SHOT, shot)
        end).id[1])
    end
end
function id!(data::hybrid_classification)
    id!(data.data_2D)
end

import Base.which
function which(feature::String, version::String=""; D::DataFrame=all_2D_features, count::Bool=false)
    if feature == "HYBRID"
        data = filter(Meta.parse("HYBRID_$(version)") => ==("YES"), D)
    elseif feature == "BASELINE"
        data = filter(Meta.parse("HYBRID_$(version)") => ==("NO"), D)
    else
        data = filter(Meta.parse(feature) => ==(1.0), D)
    end
    if !count
        return data
    else
        cnt = Dict(features_2D .=> [+(data[!, feat]...) for feat in features_2D])
        return data, cnt
    end
end
function which(features::Vector{String}, version::String=""; D::DataFrame=all_2D_features, count::Bool=false)
    data = D
    for feature in features
        data = which(feature, version; D=data)
    end
    if !count
        return data
    else
        cnt = Dict(features_2D .=> [+(data[!, feat]...) for feat in features_2D])
        return data, cnt
    end
end
function which(toks::Vector{String}, features::Union{String, Vector{String}}, version::String=""; D::DataFrame=all_2D_features, count::Bool=false) 
    D = (@subset D @byrow in(:TOK, toks))
    return which(features, version, D=D, count=count)
end
function which(tok::String, features::Union{String, Vector{String}}, version::String=""; D::DataFrame=all_2D_features, count::Bool=false) 
    D = (@subset D @byrow :TOK == tok)
    return which(features, version, D=D, count=count)
end

import Base.+
function +(variable_profiles::Vector{variable_profile}, new_feature::String)

    PROFS_ = deepcopy(variable_profiles)
    ℓ = length(variable_profiles)

    T = Vector{Vector{Float64}}(undef, ℓ)
    ρ_all = Vector{Vector{Float64}}(undef, ℓ)
    gbdf = Vector{GroupedDataFrame{DataFrame}}(undef, ℓ)

    for (k, prof) in enumerate(PROFS_)
        rename!(prof.y, ["t$(Int64(round(i * 1e6, digits=0)))" for i in prof.t])
        insertcols!(prof.y, 1, :ρ => prof.ρ) 

        T[k] = prof.t
        ρ_all[k] = prof.ρ  

        gbdf[k] = groupby(prof.y, :ρ) 
    end

    t_all = sort(unique(vcat(T...)))
    actionable_t = ["t$(Int64(round(t * 1e6, digits=0)))" for t in t_all]
    ρ_all = sort(unique(vcat(ρ_all...)))

    final = DataFrame()
    for ρi in ρ_all
        tmp = Vector{DataFrame}()
        for k in 1:ℓ
            if in(ρi, PROFS_[k].ρ)
                D = gbdf[k][(ρ=ρi,)]
                headers = names(D)
                push!(tmp, DataFrame(headers .=> mean.(eachcol(D))))
            end
        end
        dict = Dict(actionable_t .=> zeros(length(t_all)))
        dict["ρ"] = ρi
        for (t, t_label) in zip(t_all, actionable_t)
            for k in 1:ℓ
                if in(t, T[k])
                    dict["$t_label"] += tmp[k][1, "$t_label"]
                end
            end
        end
        append!(final, DataFrame(dict))
    end
    y = sort(select(final, actionable_t))
    rename!(y, actionable_t .=> ["x$i" for i in 1:length(t_all)])
    
    variable_profile(new_feature, sort(t_all), ρ_all, y, length(ρ_all), Vector{Tuple}(), DataFrame())
end

function tok_shots(feature::Union{String, Vector{String}}, version::String="")
    Tuple.(eachrow(which(feature, version)[!, [:tok, :shot]]))
end
function tok_shots(D::Union{DataFrame, SubDataFrame}) 
    Tuple.(eachrow(D[!, [:tok, :shot]]))
end

function find_2D(id::String; metadata=false)
    tok, shot = (@subset id_codes_2D @byrow :id == id)[1, [:TOK, :SHOT]]
    if !metadata
        profiles(tok, shot)
    elseif metadata
        DataFrame(:id => id, :TOK => tok, :SHOT => shot, :HYBRID => "")
    end
end

function find_0D(id::String; metadata=false) 
    D = original_space(DATA_) 
    if !metadata
        return @subset D @byrow :id == "$id"
    elseif metadata
        (@subset D @byrow :id == "$id")[!, [:id, :TOK, :SHOT, :HYBRID]]
    end
end

function classification_v1!(data_0D::data0D)
    df = data_0D.data.original_space
    class = df[!, [:id, :TOK, :SHOT, :HYBRID]]
    class.class_ID .= "original classification"
    class.comments .= ""
    data_0D.hybrid = @subset class @byrow :HYBRID == "YES"
    data_0D.baseline = @subset class @byrow :HYBRID == "NO"
    data_0D.tbd = @subset class @byrow :HYBRID == "UNKNOWN"
    return data_0D 
end
function classification_v1!(data_2D::data2D; permanant::DataFrame=all_2D_features)
    class = DataFrame()
    for (id, tok, shot, ITPA_data) in eachrow(id_codes_2D)
        if in(shot, vcat(AUG_baseline, AUG_ITER_baseline, JET_baseline))
            append!(class, DataFrame(:id => id, 
                                    :TOK => tok, 
                                    :SHOT => shot, 
                                    :HYBRID => "NO"))
            continue
        elseif in(shot, vcat(AUG_hybrid, JET_hybrid))
            append!(class, DataFrame(:id => id, 
                                    :TOK => tok, 
                                    :SHOT => shot, 
                                    :HYBRID => "YES"))
            continue
        elseif ITPA_data == 1
            D = find_0D(id, metadata=true)
            append!(class, DataFrame(D[1, :]))
            continue
        else
            append!(class, DataFrame(:id => id, 
                                    :TOK => tok, 
                                    :SHOT => shot, 
                                    :HYBRID => "UNKNOWN"))
            continue
        end
    end
    class.class_ID .= "original classification"
    class.comments .= ""

    data_2D.hybrid = @subset class @byrow :HYBRID == "YES"
    data_2D.baseline = @subset class @byrow :HYBRID == "NO"
    data_2D.tbd = @subset class @byrow :HYBRID == "UNKNOWN"
    
    permanant.HYBRID_v1 .= "UNKNOWN"
    @eachrow! permanant :HYBRID_v1 = in(:id, data_2D.hybrid.id) ? "YES" : :HYBRID_v1 
    @eachrow! permanant :HYBRID_v1 = in(:id, data_2D.baseline.id) ? "NO" : :HYBRID_v1 
end
function classification_v1!()
    classification_v1!(DATA_.data_0D)
    classification_v1!(DATA_.data_2D)
    DATA_
end

function classification_v1!(data_0D::data0D)
    df = data_0D.data.original_space
    class = df[!, [:id, :TOK, :SHOT, :HYBRID]]
    class.class_ID .= "original classification"
    class.comments .= ""
    data_0D.hybrid = @subset class @byrow :HYBRID == "YES"
    data_0D.baseline = @subset class @byrow :HYBRID == "NO"
    data_0D.tbd = @subset class @byrow :HYBRID == "UNKNOWN"
    return data_0D 
end
function classification_v1!(data_2D::data2D; permanant::DataFrame=all_2D_features)
    class = DataFrame()
    for (id, tok, shot, ITPA_data) in eachrow(id_codes_2D)
        if in(shot, vcat(AUG_baseline, AUG_ITER_baseline, JET_baseline))
            append!(class, DataFrame(:id => id, 
                                    :TOK => tok, 
                                    :SHOT => shot, 
                                    :HYBRID => "NO"))
            continue
        elseif in(shot, vcat(AUG_hybrid, JET_hybrid))
            append!(class, DataFrame(:id => id, 
                                    :TOK => tok, 
                                    :SHOT => shot, 
                                    :HYBRID => "YES"))
            continue
        elseif ITPA_data == 1
            D = find_0D(id, metadata=true)
            append!(class, DataFrame(D[1, :]))
            continue
        else
            append!(class, DataFrame(:id => id, 
                                    :TOK => tok, 
                                    :SHOT => shot, 
                                    :HYBRID => "UNKNOWN"))
            continue
        end
    end
    class.class_ID .= "original classification"
    class.comments .= ""

    data_2D.hybrid = @subset class @byrow :HYBRID == "YES"
    data_2D.baseline = @subset class @byrow :HYBRID == "NO"
    data_2D.tbd = @subset class @byrow :HYBRID == "UNKNOWN"
    
    permanant.HYBRID_v1 .= "UNKNOWN"
    @eachrow! permanant :HYBRID_v1 = in(:id, data_2D.hybrid.id) ? "YES" : :HYBRID_v1 
    @eachrow! permanant :HYBRID_v1 = in(:id, data_2D.baseline.id) ? "NO" : :HYBRID_v1 
end
function classification_v1!()
    classification_v1!(DATA_.data_0D)
    classification_v1!(DATA_.data_2D)
    DATA_
end

import MLJ.training_parition
function training_partition(labelled_data::OrderedDict, labels::Vector{String}, k::Int=2; S::Int=123)
    ℓ = countmap(labelled_data |> values)

    i = 0
    TD_ind = Vector{Int}(undef, 0)
    for label in labels
        TD_ind = vcat(TD_ind, sample(Random.seed!(S), i+1:i+ℓ[label], k, replace=false))
        i += ℓ[label]
    end
    
    return [(in(i, TD_ind) ? true : false) for i in 1:length(labelled_data)]
end
function training_partition(training::DTW_results, k::Int=2; S::Int=123)
    y = training.metadata.y
    ℓ = countmap(y)

    i = 0
    TD_ind = Vector{Int}(undef, 0)
    for yi in y
        TD_ind = vcat(TD_ind, sample(Random.seed!(S), i+1:i+ℓ[yi], k, replace=false))
        i += ℓ[yi]
    end
    
    return [(in(i, TD_ind) ? true : false) for i in 1:length(y)]
end
function training_partition(labelled_data::OrderedDict, labels::Vector{String}, k::Vector{Int}=[2]; S::Int=123)
    ℓ = countmap(labelled_data |> values)
    @assert length(labels) == length(k)

    i = 0
    TD_ind = Vector{Int}(undef, 0)
    for (label, ki) in zip(labels, k)
        TD_ind = vcat(TD_ind, sample(Random.seed!(S), i+1:i+ℓ[label], ki, replace=false))
        i += ℓ[label]
    end
    
    return [(in(i, TD_ind) ? true : false) for i in 1:length(labelled_data)]
end
function training_partition(training::DTW_results, k::Dict; S::Int=123)
    meta = deepcopy(training.metadata)
    meta.ind = 1:size(meta, 1)
    ℓ = countmap(meta.y)
    @assert length(ℓ) == length(k)
    @assert sort(collect(keys(ℓ))) == sort(collect(keys(k)))

    TD_ind = Vector{Int}(undef, 0)
    for ((label, ki)) in k
        df = @subset meta @byrow :y == label
        TD_ind = vcat(TD_ind, sample(Random.seed!(S), df.ind, ki, replace=false))
    end
    
    return [(in(i, TD_ind) ? true : false) for i in 1:length(meta.ind)]
end
function training_partition(training::DataFrame, fraction::Float64; S::Int=123)
    @assert fraction <= 1.0 "fraction needs to be a number between 0 and 1"
    
    meta = deepcopy(training)
    meta.ind = 1:size(meta, 1)
    
    ℓ = countmap(meta[!, :label])
    labels = unique(String.(meta[!, :label]))

    TD_ind = Vector{Int}(undef, 0)
    for label in labels
        df = @subset meta @byrow :label == label

        ki = Int(ceil(ℓ[label]*fraction))
        TD_ind = vcat(TD_ind, sample(Random.seed!(S), df.ind, ki, replace=false))
    end
    
    return [(in(i, TD_ind) ? true : false) for i in 1:length(meta.ind)]
end

function hyper_parameter_search(data::DTW_hyp, labelled_data::OrderedDict{Tuple{String, Int64}, String}, k::Union{Vector{Int}, Int}; interesting::String="", N::Int=5, metric::String, max::Bool=true, count::Int=7, cost::Bool=false)
    if cost
        return hyper_parameter_search_C2(data, labelled_data, k; interesting=interesting, N=N, metric=metric, max=max, count=count)
    end
    labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

    cos_arr = Array(data.cosine_cost[:, 2:end])
    ft_arr = Array(data.flat_top_cost[:, 2:end])

	shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
    labels = collect(labelled_data |> values) |> unique

    cnt, ACC, δ_acc, cos_med, ft_med = 1, 0., 100, 0., 0.

    while (cnt < count) && (ACC !== 1.) && (δ_acc > 0.0001)
        if cnt == 1
            Cs = collect(Iterators.product( 
                quantile(cos_arr[labelled_ind, labelled_ind][:], [0.05, 0.15, 0.35, 0.5, 0.65, 0.85, 0.95]),
                quantile(ft_arr[labelled_ind, labelled_ind][:], [0.05, 0.15, 0.35, 0.5, 0.65, 0.85, 0.95])
            ))
        else
            Cs = collect(Iterators.product(
                range(cos_med/2, cos_med+(cos_med/2), length=3+(2*cnt)),
                range(ft_med/2, ft_med+(ft_med/2), length=3+(2*cnt))
            ))
        end

        nc, nft = size(Cs)
        hyp_search = zeros(nc, nft, 3)
        acc_int = [zeros(nc, nft) for _ in 1:nthreads()+1]

        Threads.@threads for S in ProgressBar(1:N)
            train_ind = training_partion(labelled_data, labels, k, S=S)
            for (i, j) in Iterators.product(1:nc, 1:nft)
                CC = Cs[i, j][1]
                FTC = Cs[i, j][2]

                X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ CC).^2) .* 
                    exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ FTC).^2)

                K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
                model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=1.0)

                KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
                ỹ, _ = svmpredict(model, KK)
                acc_int[threadid()][i, j] += emp(metric)(ỹ, labelled_y[Not(train_ind)])
                hyp_search[i, j, 2] = CC
                hyp_search[i, j, 3] = FTC
            end
        end
        hyp_search[:, :, 1] = map(+, acc_int...)
        hyp_search[:, :, 1] ./= N
    
        if max
            ext = maximum(hyp_search[:, :, 1])
            println("max = ", ext)
        else
            ext = minimum(hyp_search[:, :, 1]) 
            println("min = ", ext)
        end
        inds = findall(i -> i == ext, hyp_search[:, :, 1])
        println(inds)
        ind = sample(Random.seed!(123), inds, 1)[1]

        ext = round(ext, digits=5)
        (cos_med, ft_med) = (hyp_search[ind, 2], hyp_search[ind, 3])
        println(ext, ": ($cos_med, $ft_med)")
        δ_acc = abs(-(ACC, ext))

        ACC = ext
        cnt += 1
    end

	X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ (cos_med)).^2) .* 
        exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ (ft_med)).^2) 

    K = X[labelled_ind, labelled_ind]
    model = svmtrain(K, labelled_y, kernel=Kernel.Precomputed, cost=1.0)
    println("# of SV's: ", model.SVs.nSV)

    KK = X[labelled_ind, Not(labelled_ind)]
    ỹ, confidence = svmpredict(model, KK)

    res = DataFrame(hcat(
            data.tok_shots[Not(labelled_ind)], ỹ
        ), [:ts, :predict]
    )
    if interesting !== ""
        show(stdout, "text/plain", (@subset res @byrow :predict == interesting))
    end
    hyperparameters = (cos_med, ft_med, 1.0)
    return hyper_parameters(data,
                        labelled_data,
                        model,
                        res,
                        confidence,
                        KK,
                        hyperparameters,
                        DataFrame(),
                        ACC)
end
function hyper_parameter_search_C(data::DTW_hyp, labelled_data::OrderedDict{Tuple{String, Int64}, String}, k::Union{Vector{Int}, Int}; interesting::String="", N::Int=5, metric::String, max::Bool=true, count::Int=7)
    labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

    cos_arr = Array(data.cosine_cost[:, 2:end])
    ft_arr = Array(data.flat_top_cost[:, 2:end])

	shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
    labels = collect(labelled_data |> values) |> unique

    cnt, ACC, δ_acc, cos_med, ft_med = 1, 0., 100, 0., 0.

    while (cnt < count) && (ACC !== 1.) && (δ_acc > 0.0001)
        if cnt == 1
            Cs = collect(Iterators.product( 
                quantile(cos_arr[labelled_ind, labelled_ind][:], [0.05, 0.15, 0.35, 0.5, 0.65, 0.85, 0.95]),
                quantile(ft_arr[labelled_ind, labelled_ind][:], [0.05, 0.15, 0.35, 0.5, 0.65, 0.85, 0.95])
            ))
        else
            Cs = collect(Iterators.product(
                range(cos_med/2, cos_med+(cos_med/2), length=3+(2*cnt)),
                # range(mag_med/2, mag_med*2, length=10),
                range(ft_med/2, ft_med+(ft_med/2), length=3+(2*cnt))
            ))
        end
        
        nc, nft = size(Cs)
        hyp_search = zeros(nc, nft, 3)
        acc_int = [zeros(nc, nft) for _ in 1:nthreads()]
        Threads.@threads for S in ProgressBar(1:N)
            # Random.seed!(S)
            train_ind = training_partion(labelled_data, labels, k, S=S)
            for (i, j) in Iterators.product(1:nc, 1:nft)
                CC = Cs[i, j][1]
                FTC = Cs[i, j][2]

                X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ CC).^2) .* 
                    exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ FTC).^2)
                    

                K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
                if !issymmetric(K)
                    break
                end
                model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=1.0)

                KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
                ỹ, _ = svmpredict(model, KK)
                acc_int[threadid()][i, j] += emp(metric)(ỹ, labelled_y[Not(train_ind)])
                hyp_search[i, j, 2] = CC
                hyp_search[i, j, 3] = FTC
            end
        end
        hyp_search[:, :, 1] = map(+, acc_int...)
        hyp_search[:, :, 1] ./= N
        # hyp_search |> display
    
        if max
            ext = maximum(hyp_search[:, :, 1])
            println("max = ", ext)
        else
            ext = minimum(hyp_search[:, :, 1]) 
            println("min = ", ext)
        end
        inds = findall(i -> i == ext, hyp_search[:, :, 1])
        println(inds)
        ind = sample(Random.seed!(123), inds, 1)[1]

        ext = round(ext, digits=5)
        # last_cos_med, last_mag_med = cos_med, mag_med
        (cos_med, ft_med) = (hyp_search[ind, 2], hyp_search[ind, 3])
        println(ext, ": ($cos_med, $ft_med)")
        δ_acc = abs(-(ACC, ext))

        ACC = ext
        # no_good_values = length(inds)
        cnt += 1
    end
    begin
        Cs = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1., 3., 5., 10., 30., 50., 100.]
        
        l = length(Cs)
        hyp_search = zeros(l, 5)
        acc_int = [zeros(l, 4) for _ in 1:nthreads()]

        Threads.@threads for S in ProgressBar(1:N)
            train_ind = training_partion(labelled_data, labels, k, S=S)
            for i in 1:l
                C = Cs[i]

                X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
                    exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ ft_med).^2)

                K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
                if !issymmetric(K)
                    break
                end
                model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C)

                KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
                ỹ, _ = svmpredict(model, KK)

                acc_int[threadid()][i, 1] += BalancedAccuracy()(ỹ, labelled_y[Not(train_ind)])
                acc_int[threadid()][i, 2] += emp(metric)(ỹ, labelled_y[Not(train_ind)])
                acc_int[threadid()][i, 3] += mean(Int.(model.SVs.nSV) ./ k)
                acc_int[threadid()][i, 4] += Int.(model.SVs.nSV)[end] / k[end]
            end
        end
        hyp_search[:, 1] = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1., 3., 5., 10., 30., 50., 100.]
        hyp_search[:, 2:5] = map(+, acc_int...)
        hyp_search[:, 2:5] ./= N
    end
    X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ (cos_med)).^2) .* 
        exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ (ft_med)).^2) 

    hyperparameters_meta = DataFrame(hyp_search, [:C, :BA, :μSVs, :μSV_star, :metric])
    C = hyperparameters_meta[findall(i -> i == minimum(hyperparameters_meta.metric), hyperparameters_meta.metric), :C][1]
    
    println("($cos_med, $ft_med, $C)") 
    K = X[labelled_ind, labelled_ind]
    if !issymmetric(K)
        return missing
    end
    model = svmtrain(K, labelled_y, kernel=Kernel.Precomputed, cost=C)
    println("# of SV's: ", model.SVs.nSV)

    KK = X[labelled_ind, Not(labelled_ind)]
    ỹ, confidence = svmpredict(model, KK)
    res = DataFrame(hcat(
            data.tok_shots[Not(labelled_ind)], ỹ
        ), [:ts, :predict]
    )
    if interesting !== ""
        show(stdout, "text/plain", (@subset res @byrow :predict == interesting))
    end

    hyperparameters = (cos_med, ft_med, C)
    return hyper_parameters(data,
                        labelled_data,
                        model,
                        res,
                        confidence,
                        KK,
                        hyperparameters,
                        hyperparameters_meta,
                        ACC)
end
function hyper_parameter_search_C2(data::DTW_hyp, labelled_data::OrderedDict{Tuple{String, Int64}, String}, k::Union{Vector{Int}, Int}; interesting::String="", N::Int=5, metric::String, max::Bool=true, count::Int=7)
    labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

    cos_arr = Array(data.cosine_cost[:, 2:end])
    ft_arr = Array(data.flat_top_cost[:, 2:end])

	shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
    labels = collect(labelled_data |> values) |> unique
    
    Cchoices = [0.01, 0.05, 0.1, 0.5, 1., 5., 10., 50., 100.]
    l = length(Cchoices)
    best_of_C = zeros(l, 7)
    best_of_acc = [zeros(length(labels)) for _ in 1:l]

    hyperparameters_meta = DataFrame()

    cnt, ACC, δ_acc, cos_med, ft_med = 1, 0., 100, 0., 0.
    for (n, C) in enumerate(Cchoices)
        cnt, ACC, δ_acc, cos_med, ft_med = 1, 0., 100, 0., 0.
        best_of_C[n, 1] = C

        while (cnt < count) && (ACC !== 1.) && (δ_acc > 0.0001)
            if cnt == 1
                Cs = collect(Iterators.product( 
                    quantile(cos_arr[labelled_ind, labelled_ind][:], [0.05, 0.15, 0.35, 0.5, 0.65, 0.85, 0.95]),
                    quantile(ft_arr[labelled_ind, labelled_ind][:], [0.05, 0.15, 0.35, 0.5, 0.65, 0.85, 0.95])
                ))
            else
                Cs = collect(Iterators.product(
                    range(cos_med/2, cos_med+(cos_med/2), length=3+(2*cnt)),
                    # range(mag_med/2, mag_med*2, length=10),
                    range(ft_med/2, ft_med+(ft_med/2), length=3+(2*cnt))
                ))
            end
            
            nc, nft = size(Cs)
            hyp_search = zeros(nc, nft, 3)

            acc_int = [zeros(nc, nft, 4) for _ in 1:nthreads()]
            acc_elem = [[zeros(length(labels)) for _ in 1:nc, _ in 1:nft] for _ in 1:nthreads()]

            Threads.@threads for S in ProgressBar(1:N)
                train_ind = training_partion(labelled_data, labels, k, S=S)
                for (i, j) in Iterators.product(1:nc, 1:nft)
                    CC = Cs[i, j][1]
                    FTC = Cs[i, j][2]

                    X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ CC).^2) .* 
                        exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ FTC).^2)

                    K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
                    model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C)

                    KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
                    ỹ, _ = svmpredict(model, KK)

                    int_D = DataFrame(:pred => ỹ, :actual => labelled_y[Not(train_ind)])
                    gbdf = groupby(int_D, :actual)
                    
                    acc_int[threadid()][i, j, 1] += emp(metric)(ỹ, labelled_y[Not(train_ind)])
                    acc_int[threadid()][i, j, 2] += BalancedAccuracy()(ỹ, labelled_y[Not(train_ind)])
                    acc_int[threadid()][i, j, 3] += mean(Int.(model.SVs.nSV) ./ k)
                    acc_int[threadid()][i, j, 4] += Int.(model.SVs.nSV)[end] / k[end]
                    
                    acc_elem[threadid()][i, j] += [Accuracy()(gbdf[(key,)].pred, gbdf[(key,)].actual) for key in labels]
                    hyp_search[i, j, 2] = CC
                    hyp_search[i, j, 3] = FTC
                end
            end
            acc_int = map(+, acc_int...)
            acc_int ./= N

            acc_elem = map(+, acc_elem...)
            acc_elem ./= N

            hyp_search[:, :, 1] = acc_int[:, :, 1]

            if max
                ext = maximum(hyp_search[:, :, 1])
                println("max = ", ext)
            else
                ext = minimum(hyp_search[:, :, 1]) 
                println("min = ", ext)
            end

            inds = findall(i -> i == ext, hyp_search[:, :, 1])
            println(inds)
            ind = sample(Random.seed!(123), inds, 1)[1]

            ext = round(ext, digits=5)
            best_of_C[n, 2:5] = acc_int[ind, :]
            best_of_C[n, 6:7] = [hyp_search[ind, 2], hyp_search[ind, 3]]
            best_of_acc[n] = acc_elem[ind]

            (cos_med, ft_med) = (hyp_search[ind, 2], hyp_search[ind, 3])
            println(ext, ": ($cos_med, $ft_med)")
            δ_acc = abs(-(ACC, ext))

            ACC = ext
            # no_good_values = length(inds)
            cnt += 1
        end
    end

    if max
        ext = maximum(best_of_C[:, 2])
        println("max = ", ext)
    else
        ext = minimum(best_of_C[:, 2]) 
        println("min = ", ext)
    end
    inds = findall(i -> i == ext, best_of_C[:, 2])
    ind = sample(Random.seed!(123), inds, 1)[1]

    hyperparameters_meta = DataFrame(best_of_C, [:C, :metric, :BA, :μSVs, :μSV_star, :cos, :ft])
    insertcols!(hyperparameters_meta, 4, :class_acc => best_of_acc)
    (cos_med, ft_med, C) = Array(hyperparameters_meta[ind, [:cos, :ft, :C]])

    X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ cos_med).^2) .*
        exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ ft_med).^2)
    
    K = X[labelled_ind, labelled_ind]
    model = svmtrain(K, labelled_y, kernel=Kernel.Precomputed, cost=C)
    println(labels)
    println("# of SV's: ", model.SVs.nSV)

    KK = X[labelled_ind, Not(labelled_ind)]
    ỹ, confidence = svmpredict(model, KK)
    res = DataFrame(hcat(
            data.tok_shots[Not(labelled_ind)], ỹ
        ), [:ts, :predict]
    )
    if interesting !== ""
        show(stdout, "text/plain", (@subset res @byrow :predict == interesting))
    end

    hyperparameters = (cos_med, ft_med, C)
    return hyper_parameters(data,
                        labelled_data,
                        model,
                        res,
                        confidence,
                        KK,
                        hyperparameters,
                        hyperparameters_meta,
                        ACC)
end

# 
function classify!(data::DTW_hyp, labelled_data::OrderedDict{Tuple{String, Int64}, String}, hyperparameters::Tuple; interesting::String="")

    labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

    X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ (hyperparameters[1])).^2) .* 
        exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ (hyperparameters[2])).^2) 

    K = X[labelled_ind, labelled_ind]
    model = svmtrain(K, labelled_y, kernel=Kernel.Precomputed, cost=hyperparameters[3])
    println("# of SV's: ", model.SVs.nSV)

    KK = X[labelled_ind, Not(labelled_ind)]
    # ỹ, _ = KNN(KK, labelled_y)
    ỹ, confidence = svmpredict(model, KK)

    res = DataFrame(hcat(
            data.tok_shots[Not(labelled_ind)], ỹ
        ), [:ts, :predict]
    )
    if interesting !== ""
        show(stdout, "text/plain", (@subset res @byrow :predict == interesting))
    end
    return res, confidence, KK, model
end

function find_important_training!(HYP::hyper_parameters, ts::Tuple{String, Int})
    HYP.kernel
    shots = [j for (i, j) in comp_CH[para].results.ts]
end