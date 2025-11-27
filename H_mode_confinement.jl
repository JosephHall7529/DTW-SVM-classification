abstract type confinement_data end
abstract type workflow end

mutable struct original_space_ <: confinement_data
    ELMy::DataFrame
    ITER_like::DataFrame
    JET_only::DataFrame 
    # hybrid::DataFrame
    function original_space_()
        new()
    end
end

mutable struct log_space_ <: confinement_data
    ELMy::DataFrame
    ITER_like::DataFrame
    JET_only::DataFrame 
    # hybrid::DataFrame
    function log_space_()
        new()
    end
end

mutable struct H_mode_data <: workflow
    csv::DataFrame
    modified::DataFrame
    reduced::DataFrame
    original_space::original_space_
    log_space::log_space_
    database::Symbol

    function dimensionless(df::DataFrame)
        data = deepcopy(df)

        function β(col1, col2)
            β = Vector{Union{Missing, Float64}}(undef, length(col1))
            for (n, (C1, C2)) in enumerate(zip(col1, col2))
                if C1 === missing
                    if C2 === missing
                        β[n] = missing
                    else
                        β[n] = C2
                    end
                else
                    if C2 === missing
                        β[n] = C1
                    else
                        β[n] = mean([C1, C2])
                    end
                end
            end
            return β
        end

        @transform! data @passmissing @byrow :EPS = :AMIN / :RGEO
        @transform! data @passmissing @byrow :LCOULOMB = 37.8 .- log.((:NEL .* 0.88).^0.5 ./ (:TAV .* 1e-3))
        @transform! data @passmissing @byrow :OMEGACYCL = abs.(:BT) ./ :MEFF
        @transform! data @passmissing @byrow :QCYL5 = 5e6 .* abs.(:BT ./ :IP) .* :RGEO .* :EPS .^ 2 * :KAPPAA
        @transform! data @passmissing @byrow :TAUBOHM = :OMEGACYCL * :TAUTH
        @transform! data @passmissing @byrow :RHOSTAR = 1.44e-4 .* (:MEFF .* :TAV)^0.5 ./ (abs.(:BT) .* :RGEO .* :EPS)
        @transform! data @passmissing @byrow :BETASTAR = 1.68e-4 .* :WTH ./ (:BT .^2 .* :VOL)
        @transform! data @passmissing @byrow :NUSTAR = 1e-17 .* (:NEL .* 0.88) .* :LCOULOMB .* (:QCYL5 .* :RGEO .^2.5) ./ (:TAV .^2 .* (:RGEO .* :EPS) .^1.5)
        @transform! data @passmissing @byrow :Hugill = (5 .* :AMIN.^2 .* (:NEL .* 1e-19))./(abs(:IP .* 1e-6))

        @transform! data @passmissing @byrow :ng = 10 .* abs.(:IP .* 1e-6) ./ (π .* :AMIN .^2)
        @transform! data @passmissing @byrow :nnorm = (:NEL .* 1e-19) ./ :ng
        @transform! data @passmissing @byrow :BETAN = (:BETASTAR .* :AMIN .* abs.(:BT)) ./ abs.(:IP .* 1e-6)

        @transform! data @byrow @passmissing :τE_98 = 0.0562 .* (abs.((:IP .* 1e-6)) .^ 0.93) .* (abs.(:BT) .^ 0.15) .* ((:NEL .* 1e-19) .^ 0.41) .* ((:PLTH .* 1e-6) .^ -0.69) .* (:RGEO .^ 1.97) .* (:KAREA .^ 0.78) .* (:EPS .^ 0.58) .* (:MEFF .^ 0.19)  
        @transform! data @byrow @passmissing :τE_NI = 0.068 .* (abs.((:IP .* 1e-6)) .^ 0.76) .* (abs.(:BT) .^ 0.32) .* ((:NEL .* 1e-19) .^ 0.44) .* ((:PLTH .* 1e-6) .^ -0.76) .* (:RGEO .^ 2.2) .* (:KAREA .^ 0.56) .* (:EPS .^ 0.79) .* (:MEFF .^ 0.13)  
        @transform! data @byrow @passmissing :τE_20 = 0.058 .* (abs.((:IP .* 1e-6)) .^ 0.98) .* (abs.(:BT) .^ 0.22) .* ((:NEL .* 1e-19) .^ 0.24) .* ((:PLTH .* 1e-6) .^ -0.67) .* (:RGEO .^ 1.71) .* ((1 .+ :DELTA) .^ 0.36) .* (:KAREA .^ 0.8) .* (:EPS .^ 0.35) .* (:MEFF .^ 0.2)  
        @transform! data @byrow @passmissing :τE_20IL = 0.067 .* (abs.((:IP .* 1e-6)) .^ 1.29) .* (abs.(:BT) .^ -0.13) .* ((:NEL .* 1e-19) .^ 0.15) .* ((:PLTH .* 1e-6) .^ -0.64) .* (:RGEO .^ 1.19) .* ((1 .+ :DELTA) .^ 0.56) .* (:KAREA .^ 0.67) .* (:MEFF .^ 0.3)  
        @transform! data @byrow @passmissing :τn98 = :TAUTH ./ :τE_98
        @transform! data @byrow @passmissing :τnNI = :TAUTH ./ :τE_NI
        @transform! data @byrow @passmissing :τn20 = :TAUTH ./ :τE_20
        @transform! data @byrow @passmissing :τn20IL = :TAUTH ./ :τE_20IL

        data.β = β(data.BEPDIA, data.BEPMHD)
        @transform! data @passmissing @byrow :βt_3p5 = (3.5 .* (:IP .* 1e-6)) ./ (:AMIN .* abs.(:BT))
        @transform! data @passmissing @byrow :βt_VA = 1.68e-4 .* (:WTH) ./ (abs.(:BT) .^2 .* :VOL)
        @transform! data @passmissing @byrow :βp_3p5 = (3.5^2 .* :KAREA)/(4 .* :βt_VA)
        # @transform! data @passmissing @byrow :βN = (:βt_VA .* :AMIN .* abs.(:BT)) ./ abs.(:IP)

        # Bootstrap, needs refinement
        @transform! data @passmissing @byrow :Ip_bs = 0.65 .* (:EPS .^ 0.5) .* :β .* (abs.(:IP) .* 1e-6)
        @transform! data @passmissing @byrow :fbs = :Ip_bs ./ abs.(:IP * 1e-6)
        
        @transform! data @passmissing @byrow :FoM = (:BETAN .* :τn98) ./ (:Q95.^2)
        @transform! data @passmissing @byrow :Fluence_FoM = (:BETAN.^(3.5) .* abs(:BT).^(3.5) .* :AMIN .^(2.5)) ./ (:fbs .^(1.5))
        return data
    end

    function data_modify(data::DataFrame)
        select!(data, :, :IP => (x -> abs.(x .* 1e-6)) => :IP, :) # converting to MegaAmps
        select!(data, :, :PLTH => (x -> x .* 1e-6) => :PLTH, :) # converting to Mega# Converting to MegaWatts
        select!(data, :, :PL => (x -> x .* 1e-6) => :PL, :) # converting to Mega# Converting to MegaWatts
        select!(data, :, :NEL => (x -> x .* 1e-19) => :NEL, :)  # converting too 10^19 m-3 
        select!(data, :, :WTH => (x -> x .* 1e-6) => :WTH, :)  # converting too MJ
        return data
    end
    
    function H_mode_data(csv::AbstractString)
        df = CSV.read(csv, DataFrame, header=true, stringtype=String)
        df = unique_data(df, unique_comparison_parameters).corrected
        df = dimensionless(df)
        new(df, data_modify(df))
    end
    function H_mode_data(df::DataFrame)
        df = unique_data(df, unique_comparison_parameters).corrected
        df = dimensionless(df)
        new(df, data_modify(df))
    end
    function H_mode_data()
        new()
    end
end

function id(H::H_mode_data)
    return H.reduced.id
end

import Base.abs
function abs(T::Tuple)
    return abs.(T)
end
function abs(V::Vector)
    return abs.(V)
end

import Base.log
function log(T::Tuple)
    log.(T)
end
function log(V::Vector)
    log.(V)
end
function log(H::H_mode_data, field::Symbol)
    df = getfield(H.original_space, field)
    df_abs = select(df, :, feature_parameters .=> ByRow(abs), :, renamecols=false)
    return select(df_abs, :, feature_parameters .=> ByRow(log), :, renamecols=false)
end

function H_mode_reduce!(H::H_mode_data, flag::Symbol)
    df = H.modified
    params = union(parameters, parameters_dimensionless)
    if flag === :none
        df = select(df, params)
    else
        df = @subset(df, $flag .== 1)
    end
    df = dropmissing(df, vcat(regression_parameters))
    H.reduced = df
end

function ELMy_reduce!(H::H_mode_data)
    df = H.reduced
    df = @subset(df, @byrow in(:PHASE, ELMy_indicators))
    H.original_space.ELMy = df
    H.log_space.ELMy = log(H, :ELMy)
end

function H_mode_fill!(H::H_mode_data, flag::Symbol)
    H_mode_reduce!(H, flag)
    original = original_space_()
    log = log_space_()
    H.original_space = original
    H.log_space = log
    ELMy_reduce!(H)
    # ITER_reduce!(H)
    # JET_reduce!(H)
    H.database = flag
end

function cardinality_metadata(H::H_mode_data, field::Symbol; name::Symbol=:cardinality)
    no_of_points = DataFrame()
    df = getfield(H.original_space, field)

    gbdf, ind = grouping(df, :TOK)

    for id in ind
        ℓ = size(gbdf[id], 1)
        df = DataFrame(:TOK => id[1], name => ℓ)
        no_of_points = vcat(no_of_points, df)
    end
    df = DataFrame(:TOK => "Total", name => sum(no_of_points[!, name]))
    no_of_points = vcat(no_of_points, df)
    return no_of_points
end

function correlation_metadata(H::H_mode_data, fields::Symbol=:ELMy, space::Symbol=:log_space_; parameters::Vector=regression_predictor, save=false, kwargs...)
    D = getfield(H, space)
    ℓ = length(parameters)

    reduced_1 = getfield(D, fields)[!, parameters] |> Array
    data_1 = StatsBase.standardize(ZScoreTransform, reduced_1, dims=1)

    C1 = cor(data_1)

    df = DataFrame(round.(C1; kwargs...), parameters)
    insertcols!(df, 1, :_ => parameters)
    return df   
end
function correlation_metadata(H::H_mode_data, fields::Vector{Symbol}=[:ELMy, :ITER_like], space::Symbol=:log_space_; parameters::Vector=regression_predictor, save=false, kwargs...)
    fields = sort(fields)
    D = getfield(H, space)
    ℓ = length(parameters)

    reduced_1 = getfield(D, fields[1])[!, parameters] |> Array
    data_1 = StatsBase.standardize(ZScoreTransform, reduced_1, dims=1)

    reduced_2 = getfield(D, fields[2])[!, parameters] |> Array
    data_2 = StatsBase.standardize(ZScoreTransform, reduced_2, dims=1)

    C1 = cor(data_1)
    C2 = cor(data_2)

    correlation = [i < j ? C1[i, j] : C2[i, j] for i in 1:ℓ, j in 1:ℓ]
    df = DataFrame(round.(correlation; kwargs...), parameters)
    insertcols!(df, 1, :_ => parameters)
    if save
        dir = joinpath(meta_data_dir[1], "correlation/$(directory_flag(H))/$(space)")
        mkpath(dir)
        CSV.write(joinpath(dir, "$(*(fields, delim=:_))_correlation_matrix.csv"), df)
    end
    return df   
end

function vif_log(D::H_mode_data, field::Symbol, parameters=regression_predictor; kwargs...)
    data = getfield(D.log_space, field)[!, parameters]
    
    variations = Float64[]
    for var in parameters
        avg = mean(data[!, var])
        push!(variations, sum((data[!, var] .- avg).^2))
    end

    return DataFrame(:coefficients => (:α .* parameters),
                    field => vif(data; kwargs...),
                    :var => round.(vif(data; kwargs...) ./ variations, digits=3))
end
function vif_original(D::H_mode_data, field::Symbol, parameters=regression_predictor; kwargs...)
    data = getfield(D.original_space, field)[!, parameters]

    variations = Float64[]
    for var in parameters
        avg = mean(data[!, var])
        push!(variations, sum((data[!, var] - avg).^2))
    end

    return DataFrame(:coefficients => (:α .* parameters),
                    field => vif(data; kwargs...),
                    :var => round.(vif(data; kwargs...) ./ variations, digits=3))
end
function vif(D::H_mode_data, field::Symbol, parameters=regression_predictor; original::Bool=false, kwargs...)

    if !original 
        return vif_log(D, field; kwargs...)
    else
        return vif_original(D, field; kwargs...)
    end
end

function condition_index_log(D::H_mode_data, field::Symbol, parameters=regression_predictor; kwargs...)
    data = getfield(D.log_space, field)[!, parameters]

    return condition_index(data; kwargs...)
end
function condition_index_original(D::H_mode_data, field::Symbol, parameters=regression_predictor; kwargs...)
    data = getfield(D.original_space, field)[!, parameters]

    return DataFrame(:coefficients => (:α .* parameters),
                    field => condition_index(data; kwargs...))
end
function condition_index(D::H_mode_data, field::Symbol, parameters=regression_predictor; original::Bool=false, kwargs...)

    if !original 
        return condition_index_log(D, field, parameters; kwargs...)
    else
        return condition_index_original(D, field, parameters; kwargs...)
    end

end

function H_mode_data(csv::Union{AbstractString, DataFrame}, flag::Symbol)
    data = H_mode_data(csv)
    H_mode_fill!(data, flag)

    insertcols!(data.original_space.ELMy, :deviation => DB5_dev_labels)
    insertcols!(data.log_space.ELMy, :deviation => DB5_dev_labels)

    return data
end