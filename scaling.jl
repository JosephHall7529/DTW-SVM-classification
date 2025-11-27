abstract type regression_method end

abstract type regression_model end

struct ols <: regression_method
    method::Symbol
    function ols()
        new(:ols)
    end
end

struct regressor_variables
    parameters::Vector
    dependent::Vector
    predictors::Vector
    function regressor_variables(parameters::Vector, dependent_ind::Int64)
        new(parameters, parameters[dependent_ind], deleteat!(parameters, dependent_ind))
    end
    function regressor_variables(dependent::Vector, predictor::Vector)
        new(vcat(dependent, predictor), dependent, predictor)
    end
end

mutable struct regressor_coefficients
    names::Vector
    values::Vector
    function regressor_coefficients(predictor::regressor_variables)
        if predictor.predictors == regression_predictor
            coeffs = regression_coefficients
            new(coeffs)
        else
            coeffs = :α .* predictor.predictors
            new(vcat(:α0, coeffs))
        end
    end
    function regressor_coefficients(coeffs::Vector)
        new(coeffs)
    end
end

struct IPB98 <: regression_model
    formula::FormulaTerm
    parameters::regressor_variables
    coefficients::regressor_coefficients
    kadomtsev::Bool
    function IPB98()
        reg_var = regressor_variables(regression_dependent, regression_predictor)
        reg_coeff = regressor_coefficients(reg_var)
        new(@formula(TAUTH ~ IP + BT + NEL + PLTH + RGEO + KAREA + EPS + MEFF), reg_var, reg_coeff, false)
    end
    function IPB98(dependent::String)
        reg_var = regressor_variables([Symbol(dependent)], regression_predictor)
        reg_coeff = regressor_coefficients(reg_var)
        new(emp("@formula($(dependent) ~ IP + BT + NEL + PLTH + RGEO + KAREA + EPS + MEFF)"), reg_var, reg_coeff, false)
    end
end
struct single_machine <: regression_model
    formula::FormulaTerm
    parameters::regressor_variables
    coefficients::regressor_coefficients
    kadomtsev::Bool
    function single_machine()
        reg_var = regressor_variables(regression_dependent, [:IP, :BT, :NEL, :PLTH])
        reg_coeff = regressor_coefficients(reg_var)
        new(@formula(TAUTH ~ IP + BT + NEL + PLTH), reg_var, reg_coeff, false)
    end
    function single_machine(dependent::String)
        reg_var = regressor_variables([Symbol(dependent)], [:IP, :BT, :NEL, :PLTH])
        reg_coeff = regressor_coefficients(reg_var)
        new(emp("@formula($(dependent) ~ IP + BT + NEL + PLTH)"), reg_var, reg_coeff, false)
    end
end
struct clemente <: regression_model
    formula::FormulaTerm
    parameters::regressor_variables
    coefficients::regressor_coefficients
    kadomtsev::Bool
    function clemente()
        reg_var = regressor_variables(regression_dependent, [:IP, :BT, :NEL, :PLTH, :RGEO])
        reg_coeff = regressor_coefficients(reg_var)
        new(@formula(TAUTH ~ IP + BT + NEL + PLTH + RGEO), reg_var, reg_coeff, false)
    end
end

mutable struct Regression
    data::DataFrame
    regression_method::Symbol
    model::regression_model
    results::Union{StatsModels.TableRegressionModel, Array}
    function Regression(data::DataFrame, regression_method::ols, model::regression_model=IPB98(); kwargs...) 
        new(data, regression_method.method, model, ols(data, model; kwargs...))
    end
    # function Regression(data::DataFrame, regression_method::wls, model::regression_model=IPB98(); kwargs...)
    #     new(data, regression_method.method, model, wls(data, model; kwargs...))
    # end		
    # function Regression(data::DataFrame, regression_method::ridge, model::regression_model=IPB98(); kwargs...)
    #     new(data, regression_method.method, model, ridge(data, model; kwargs...))
    # end	
    # function Regression(data::DataFrame, regression_method::pcr, model::regression_model=IPB98(3); kwargs...)
    #     new(data, regression_method.method, model, pcr(data, model; kwargs...))
    # end	
end

function ols(data::DataFrame, model::regression_model; transform::Symbol=:none)
    df = select(data, model.parameters.parameters)
    ℓ = size(df, 1)
    if ℓ < 10
        df = repeat(df, 20)
    end
    if transform == :none
        nothing
    elseif transform == :centre
        select!(df, model.parameters.parameters .=> i -> i .- mean(i), renamecols=false)
    elseif transform == :standardize
        df = DataFrame(StatsBase.standardize(ZScoreTransform, df |> Array, dims=1), model.parameters.parameters)
    end
    result = lm(model.formula, df)
    model.coefficients.values = result |> StatsBase.coef
    model.coefficients.values 
end

function predict(D::DataFrame, results::StatsModels.TableRegressionModel, model::regression_model)
    df = Array(D[!, model.parameters.predictors])
	return [dot(StatsBase.coef(results), vcat(1.0, pnt)) for pnt in eachrow(df)]
end
function predict(D::DataFrame, results::Vector, model::regression_model)
    df = Array(D[!, model.parameters.predictors])
	return [dot(results, vcat(1.0, pnt)) for pnt in eachrow(df)]
end
function predict(D::DataFrame, results::DataFrame, model::regression_model)
    df = Array(D[!, model.parameters.predictors])
	return [dot(Array(results[!, model.coefficients.names]), vcat(1.0, pnt)) for pnt in eachrow(df)]
end