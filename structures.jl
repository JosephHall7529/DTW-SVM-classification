mutable struct variable_profile
    feature::String
    t::Vector{Float64}
    ρ::Vector{Float64}
    y::DataFrame
    ℓ::Int64
    regions::Vector{Tuple}
    stationary_times::DataFrame
end

struct stationary_variable_profile
    feature::String
    t::Vector{Float64}
    ρ::Vector{Float64}
    y::DataFrame
    ℓ::Int64
    regions::Vector{Tuple}
    stationary_times::DataFrame
end

import Base.vcat
function vcat(profile::variable_profile)
    @assert size(profile.stationary_times, 1) !== 0 "Need to run stationary!(P)"
    @assert length(unique(profile.stationary_times.stationary)) == 2 "no stationary intervals"

    st = @subset profile.stationary_times @byrow :stationary == 1
    ind = findall(in(st.t), profile.t)
    return stationary_variable_profile(profile.feature,
                                        st.t,
                                        profile.ρ,
                                        (@view profile.y[!, ind]),
                                        profile.ℓ,
                                        profile.regions,
                                        (@view profile.stationary_times[ind, :]))
end

mutable struct ITPA_profile_extract
    features::Vector{String}
    data::OrderedDict{String, variable_profile}
    tok::String
    shot::Int64
    id::String

    function ITPA_profile_extract() 
        new()
    end

    function ITPA_profile_extract(str::String)
        file=open(str, "r")

        features = Vector{String}()
        profile = OrderedDict{String, variable_profile}()
    
        while !eof(file)
            io = readuntil(file, "(Y) LABEL-")
            io = readuntil(file, ";")
            dep = String.(split(io)[1])
            if dep == "BPOL"
                dep = "BETAPOL"
            end
            push!(features, dep)
    
            io = readuntil(file, "FLAG")
            io = readuntil(file, ";")
            no_ρ = parse(Int64, split(io)[1])
    
            io = readuntil(file, "PTS-")
            io = readuntil(file, ";")
            no_t = parse(Int64, split(io)[1])
    
            io = readuntil(file, ":")
            io = readuntil(file, ";")
    
            data = [parse(Float64, y.match) for y in eachmatch(r".*?.E.{1,3}", io)]
            ρ = data[1:no_ρ]
            t = data[no_ρ+1:no_ρ+no_t] 
            y = data[no_ρ+no_t+1:end]

            arr = Float64[] 
            for (i, ti) in enumerate(t)
                append!(arr, y[no_ρ*(i-1)+1:i*no_ρ])
            end 
            D = DataFrame(reshape(arr, no_ρ, :), :auto)
            
            profile[dep] = variable_profile(dep, t, ρ, D, no_ρ, Vector{Tuple}(), DataFrame())
    
            io = readuntil(file, "*\n") 
            io = readuntil(file, "*\n") 
        end
        new(features, profile)
    end
end

mutable struct SAL_profile_extract
    diagnostics::Vector{String}
    features::Dict{String, Vector{String}}
    data::OrderedDict{String, OrderedDict{String, variable_profile}}
    shot::Int
    id::String

    function SAL_profile_extract() 
        new()
    end

    function SAL_profile_extract(shot::Int)
        features = Dict{String, Vector{String}}()
        available_diagnostics = String[]
        vp = OrderedDict{String, OrderedDict{String, variable_profile}}()
        for diagnostic in jet_diag
            println("jet: ", shot, " ", diagnostic)
            dir = files(SAL_profile_data*"/$shot/$diagnostic/", false)
            if !isempty(dir)
                push!(available_diagnostics, diagnostic)
                feat_dict = OrderedDict()
                features[diagnostic] = String[]
                for file in dir
                    feature = String(split(file, ".")[1])
                    if feature == "xlim"
                        continue
                    end
                    # println("$shot/$diagnostic/$feature")
                    push!(features[diagnostic], SAL_features_dict[feature])
                    D = CSV.read(SAL_profile_data*"/$shot/$diagnostic/$feature.csv", DataFrame)
                    if size(D, 2) == 2
                        y = DataFrame(:y => D[!, 2])
                        feat_dict[SAL_features_dict[feature]] = variable_profile(
                                                                    SAL_features_dict[feature],
                                                                    D.t,
                                                                    Vector{Float64}(),
                                                                    y,
                                                                    0,
                                                                    Vector{Tuple}(),
                                                                    DataFrame()
                        )
                    # else
                    #     y = D[!, 3:end]
                    #     ρ = emp.(names(y))
                    #     y = DataFrame(Array(y)', :auto)
                    #     ℓ = length(ρ)
                    #     # rename!(y, ["x$i" for i in 1:ℓ])
                    #     feat_dict[SAL_features_dict[feature]] = variable_profile(
                    #                                                 SAL_features_dict[feature],
                    #                                                 D.time,
                    #                                                 ρ,
                    #                                                 y,
                    #                                                 ℓ,
                    #                                                 Vector{Tuple}(),
                    #                                                 DataFrame()
                    #     )
                    end
                end
            else
                continue
            end
            vp[diagnostic] = feat_dict
        end
        new(available_diagnostics, features, vp, shot)
    end
end

mutable struct AUG_profile_extract
    diagnostics::Vector{String}
    features::Dict{String, Vector{String}}
    data::OrderedDict{String, OrderedDict{String, variable_profile}}
    shot::Int
    id::String

    function AUG_profile_extract() 
        new()
    end

    function AUG_profile_extract(shot::Int)
        features = Dict{String, Vector{String}}()
        available_diagnostics = String[]
        vp = OrderedDict{String, OrderedDict{String, variable_profile}}()
        for diagnostic in aug_diag
            dir = files(AUG_profile_data*"/$shot/$diagnostic/", false)
            if !isempty(dir)
                push!(available_diagnostics, diagnostic)
                feat_dict = OrderedDict()
                features[diagnostic] = String[]
                for file in dir
                    feature = String(split(file, ".")[1])
                    # if feature == "xlim"
                    #     continue
                    # end
                    println("$shot/$diagnostic/$feature")
                    push!(features[diagnostic], feature)
                    D = CSV.read(AUG_profile_data*"/$shot/$diagnostic/$feature.csv", DataFrame)
                    if size(D, 2) == 2
                        y = DataFrame(:y => D[!, 2])
                        feat_dict[feature] = variable_profile(feature,
                                                            D.t,
                                                            Vector{Float64}(),
                                                            y,
                                                            0,
                                                            Vector{Tuple}(),
                                                            DataFrame()
                        )
                    else
                        # y = D[!, 2:end]
                        # ρ = emp.(names(y))
                        # y = DataFrame(Array(y)', :auto)
                        # ℓ = length(ρ)
                        # # rename!(y, ["x$i" for i in 1:ℓ])
                        # feat_dict[SAL_features_dict[feature]] = variable_profile(
                        #                                             SAL_features_dict[feature],
                        #                                             D.time,
                        #                                             ρ,
                        #                                             y,
                        #                                             ℓ,
                        #                                             Vector{Tuple}(),
                        #                                             DataFrame()
                        # )
                    end
                end
            else
                continue
            end
            vp[diagnostic] = feat_dict
        end
        new(available_diagnostics, features, vp, shot)
    end
end

function ITPA_profile_extract(TOK::String, shot::Int64)
    str = "$(ITPA_profile_data)/$(TOK)/$(shot)/pr08_$(TOK)_$(shot)_2d.dat"
    @assert isfile(str)

    PE = ITPA_profile_extract(str)
    PE.tok = Unicode.normalize("$TOK", casefold=true)
    PE.shot = shot

    return PE
end

mutable struct machine_data
    tok::Union{String, Symbol}
    shots::Vector{Int64}
    profiles::Dict{Int64, Union{ITPA_profile_extract, SAL_profile_extract, AUG_profile_extract}}

    function machine_data()
        new()
    end
    function machine_data(str::String)
        str = Unicode.normalize("$str", casefold=true)

        shots = deepcopy(emp("pr08_$(str)_shots"))
        data = machine_data()
        data.tok = str
        data.shots = shots

        profiles = Dict()
        for shot in shots
            profiles[shot] = ITPA_profile_extract(str, shot)
        end

        if str == "jet"
            append!(data.shots, JET_shots)
            for shot in JET_shots
                profiles[shot] = SAL_profile_extract(shot)
            end
        end

        if str == "aug"
            append!(data.shots, AUG_shots)
            for shot in AUG_shots
                profiles[shot] = AUG_profile_extract(shot)
            end 
        end

        data.profiles = profiles
        return data
    end
end

struct global_data
    original_space::DataFrame
    log_space::DataFrame
end 

mutable struct data0D
    data::global_data
    baseline::DataFrame
    hybrid::DataFrame
    tbd::DataFrame
    function data0D()
        new()
    end
    function data0D(D::DataFrame)   
        df = deepcopy(D)
        LD = log(df)
        OD = select(D, [:id, :TOK, :SHOT], :)
        GD = global_data(OD, LD)
        new(GD, [DataFrame() for _ in 1:3]...)
    end
end

mutable struct data2D
    data::Dict{Union{String, Symbol}, machine_data}
    baseline::DataFrame
    hybrid::DataFrame
    tbd::DataFrame
    function data2D()
        data = Dict()
        for tok in pr08_toks
            data[tok] = machine_data(tok)
        end
        new(data, [DataFrame() for _ in 1:3]...)
    end
end

struct hybrid_classification
    data_0D::data0D
    data_2D::data2D
    function hybrid_classification(D::DataFrame)
        new(data0D(D), data2D())
    end
end



mutable struct DTW_hyp
    tok_shots::Vector{Tuple{String, Int}}
    features::Vector{String}
    shot_dict::Dict{Tuple, Int}
    profile_data::Dict{Tuple, Array}
    flat_top_data::Dict{Tuple, Array}
    cosine_cost::DataFrame
    magnitude_cost::DataFrame
    flat_top_cost::DataFrame

    function DTW_hyp()
        new()
    end
    function DTW_hyp(tok_shotz::Vector{Tuple{String, Int}}, features::Vector{String}; L=50, transportcost=1.1, short_IP::Bool=false) 
        
        profile_data = Dict{Tuple, Array{Float64}}()
        flat_top_data = Dict{Tuple, Array{Float64}}()

        features = sort(features)
        ℓf = length(features)
        ℓT = length(tok_shotz)
        
        dict_shot = Dict(1:ℓT .=> tok_shotz)
        shot_dict = Dict(tok_shotz .=> 1:ℓT)

        LSD = least_data_length(tok_shotz, features; short_IP=short_IP)
        for (i, ts) in enumerate(tok_shotz)
            mat_ip = Float64[]
            mat_nbi = Float64[]

            IP_step_constraint = LSD[ts]["IP"] 
            if IP_step_constraint.length > L
                IP_step_constraint = merge(IP_step_constraint, (length=L,))
            end
            NBI_step_constraint = LSD[ts]["NBI"] 
            if NBI_step_constraint.length > L
                NBI_step_constraint = merge(NBI_step_constraint, (length=L,))
            end

            constraint_feature_IP = profiles(ts..., IP_step_constraint.feat).t
            if short_IP
                constraint_time_IP = constraint_feature_IP[FTIPs[ts][1] .< constraint_feature_IP .< FTIPs[ts][2]]
            else
                constraint_time_IP = constraint_feature_IP[FTIP[ts][1] .< constraint_feature_IP .< FTIP[ts][2]]
            end
            println(ts, ": ", length(constraint_time_IP), " ", IP_step_constraint.length)
            constraint_step_IP = div(length(constraint_time_IP), IP_step_constraint.length)
            constraint_time_IP = constraint_time_IP[1:constraint_step_IP:end] 

            constraint_feature_NBI = profiles(ts..., NBI_step_constraint[1]).t
            constraint_time_NBI = constraint_feature_NBI[FTNBI[ts][1] .< constraint_feature_NBI .< FTNBI[ts][2]]
            constraint_step_NBI = div(length(constraint_time_NBI), NBI_step_constraint.length)
            constraint_time_NBI = constraint_time_NBI[1:constraint_step_NBI:end] 
            for (n, feat) in enumerate(features)
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

        cosine_cost = zeros(ℓT, ℓT)
        flat_top_cost = zeros(ℓT, ℓT)

        path_cos_ℓ = zeros(ℓT, ℓT) 
        path_ft_ℓ = zeros(ℓT, ℓT) 
        for i in ProgressBar(1:ℓT)
            for j in 1:ℓT
                data_j = profile_data[dict_shot[j]]
                data_i = profile_data[dict_shot[i]]

                COST_cos, path_cos_a, path_cos_b = dtw(data_j, data_i, CosineDist(); transportcost=transportcost)
                cosine_cost[i, j] = COST_cos
                path_cos_ℓ[i, j] = length(path_cos_a)

                data_j = flat_top_data[dict_shot[j]]
                data_i = flat_top_data[dict_shot[i]]

                COST_ft, path_ft_a, path_ft_b = dtw(data_j, data_i, Euclidean(); transportcost=transportcost)
                flat_top_cost[i, j] = COST_ft
                path_ft_ℓ[i, j] = length(path_ft_a)
            end
        end
        
        cosine_cost = round.(cosine_cost ./ path_cos_ℓ, sigdigits=5)
        flat_top_cost = round.(flat_top_cost ./ path_ft_ℓ, sigdigits=5) 
        
        df_names = df_ts_naming(tok_shotz)

        cosine_cost = DataFrame(cosine_cost[:, :], df_names)
        flat_top_cost = DataFrame(flat_top_cost[:, :], df_names) 
        insertcols!(cosine_cost, 1, :shots => tok_shotz) 
        insertcols!(flat_top_cost, 1, :shots => tok_shotz)

        new(tok_shotz, 
            features, 
            shot_dict, 
            profile_data, 
            flat_top_data,
            cosine_cost,
            DataFrame(), 
            flat_top_cost
        )
    end
end

mutable struct hyper_parameters
    data
    labelled_data::Union{Dict, OrderedDict}
    model
    results::DataFrame
    confidence::Array
    kernel::Array
    hyperparameters::Tuple
    hyperparameters_meta::DataFrame
    acc::Number
    function hyper_parameters()
        new()
    end
    function hyper_parameters(data, lts, model, res, conf, K, hyp, hyp_meta, acc)
        new(data, lts, model, res, conf, K, hyp, hyp_meta, acc)
    end
end

function hyperparameters_meta!(HYP::hyper_parameters; N::Int=50)
    HYP.hyperparameters_meta = DataFrame()
    begin
        cntmap = countmap(values(HYP.labelled_data))
        CLn, CEn, COn = cntmap["LH"], cntmap["EH"], cntmap["CO"]
        labelled_ts = collect(HYP.labelled_data |> keys)
        labelled_y = collect(HYP.labelled_data |> values)
        shot_dict = Dict([a => n for (n, a) in enumerate(HYP.data.tok_shots)])
        labelled_ind = [shot_dict[k] for k in labelled_ts]
        labels = labelled_y |> unique

        cos_, ft_, C_ = HYP.hyperparameters
        X = exp.(-(Array(HYP.data.cosine_cost[!, 2:end]) ./ (cos_)).^2) .* 
            exp.(-(Array(HYP.data.flat_top_cost[!, 2:end]) ./ (ft_)).^2) 
        
        k = [Int(div(*(CLn, 7), 10)), Int(div(*(CEn, 7), 10)), Int(div(*(COn, 7),10))]
    end

    Results = DataFrame()
    for S in 0:N
        begin
            train_ind = training_partition(HYP.labelled_data, labels, k, S=S)

            K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
            model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C_)

            KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
            ỹ, _ = svmpredict(model, KK)
            # CO_EH_ind = findall(i -> in(i, ["CO", "C-EH"]), ỹ)
            
            acc = Accuracy()(ỹ, labelled_y[Not(train_ind)])
            b_acc = BalancedAccuracy()(ỹ, labelled_y[Not(train_ind)])

            fpr = MulticlassFalsePositiveRate(;levels=["LH", "EH", "CO"], perm=[1,2,3])(ỹ, labelled_y[Not(train_ind)])
            println("accuracy = $(acc), balanced accuracy = $(b_acc), false positive rate = $(fpr)")

            SVL, SVE, SVO = model.SVs.nSV
            println(SVL, ", $(SVE), $(SVO)")
            cnt_map = countmap(ỹ)
            no_CO = haskey(cnt_map, "CO") ? cnt_map["CO"] : 0
            no_EH = haskey(cnt_map, "EH") ? cnt_map["EH"] : 0
            no_LH = haskey(cnt_map, "LH") ? cnt_map["LH"] : 0

            int = OrderedDict("rng" => S, "Accuracy" => round(acc, digits=3), "Balanced Accuracy" => round(b_acc, digits=3), "False Positive Rate" => round(fpr, digits=4),
                "#SV LH" => SVL, "#SV EH" => SVE, "#SV CO" => SVO,
                "#LH" => no_LH, "#EH" => no_EH, "#CO" => no_CO)
            append!(Results, DataFrame(int))
        end
    end
    HYP.hyperparameters_meta = Results
end 