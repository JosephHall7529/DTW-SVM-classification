using Pkg
Pkg.activate("/Users/juur_hall/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/src/")

using 
    CSV,
    DataFrames,
    OrderedCollections,
    DataFramesMeta,
    LinearAlgebra,
    MLJ,
    GLMakie,
    GLM,    
    Unicode,
    ColorSchemes,
    Random,
    LaTeXStrings,
    Combinatorics,
    DynamicAxisWarping, 
    Distances, 
    StatsBase, 
    ProgressBars,
    LIBSVM,
    ThreadTools,
    StatisticalMeasures

import DataFrames: select

colors = ColorSchemes.colorschemes[:Set2_8];
colors_15 = get(ColorSchemes.Set2_8, range(0.0, 1.0, length=18));
param_colors = vcat(get(ColorSchemes.dracula, range(0, 1, length=8))[1:4], shuffle(MersenneTwister(4), colors_15)[5:8]);

include("pre_defined_variables.jl")
include("extending_functions.jl")
include("preprocessing.jl")
include("H_mode_confinement.jl")
include("structures.jl")
include("MLJ_implementation.jl")
include("structure_functions.jl")
include("scaling.jl")
# include("stationary_2D.jl")
include("secondry.jl")
include("plotting.jl")

global_database = CSV.read("/Users/juur_hall/.data/Multi_Machine_Fusion_Data_ITER_21_03_01/DB52P3_ed_!S_good_ids.csv", DataFrame, stringtype=String)
SELDB5 = H_mode_data("/Users/juur_hall/.data/Multi_Machine_Fusion_Data_ITER_21_03_01/DB52P3_ed_!S_good_ids.csv", :SELDB5)

DATA_ = hybrid_classification(SELDB5.original_space.ELMy)
id!(DATA_)

FTNBI = flat_top_NBI(tok_shots(which(["PNBI"])), 0.6)
FTIP = flat_top_IP(tok_shots(which(["IP"])), 0.9)
FTIPs = Dict()
for (ts, time_window) in FTIP
    if ts[1] == "aug"
        FTIPs[ts] = (time_window[1], FTNBI[ts][1])
    end
end

# data_dtw = let
#     dict = OrderedDict()
#     # collect(powerset(["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "LI", "NGW"], 1, 1))
#     for para in [["IP", "PNBI"]]
#         println(para)
#         tss = tok_shots((@subset which(para) @byrow in(:tok, ["aug"])))
#         dict[naming(para)] = DTW_hyp(tss, para, L=100) 
#     end
#     dict
# end

begin
    tss = tok_shots((@subset which(["IP", "PNBI"]) @byrow in(:tok, ["aug"])))
    shots = [ts[2] for ts in tss]
    useful_2D_features = @subset all_2D_features @byrow begin
        :tok == "aug"
        in(:shot, shots)
    end
    allLH = sort(tok_shots((@subset useful_2D_features @byrow in(:current_heating, ["C-LH"]))))
    allEH = sort(tok_shots((@subset useful_2D_features @byrow in(:current_heating, ["C-EH"]))))
    allCO = sort(tok_shots((@subset useful_2D_features @byrow in(:current_heating, ["CO"]))))

    trials_D = let 
        df = DataFrame(:shots => allCO, :label => "CO")
        append!(df, DataFrame(:shots => allEH, :label => "EH"))
        append!(df, DataFrame(:shots => allLH, :label => "LH"))

        for (T, (S, CLn, CEn, COn)) in trial_dataset_meta
            println(T, " : ", "[$S, $(CLn), $(CEn), $(COn)]")
            CLchoices = sort(sample(Random.seed!(S), allLH, CLn, replace=false))
            CEchoices = sort(sample(Random.seed!(S), allEH, CEn, replace=false))
            COchoices = sort(sample(Random.seed!(S), allCO, COn, replace=false))
            println("$T: CL = [$(length(CLchoices))] \n CE = [$(length(CEchoices))] \n CO = [$(length(COchoices))] \n")

            int = DataFrame()
            for (choice, bulk) in zip([COchoices, CEchoices, CLchoices], [allCO, allEH, allLH])

                trial_ind = zeros(Int, length(bulk))
                pointer = [findall(i -> i == C, bulk)[1] for C in choice]
                trial_ind[pointer] .= 1

                dict = OrderedDict(:shots => [ts for ts in bulk], T => trial_ind)
                D = DataFrame(dict)
                append!(int, D)
            end
            df = innerjoin(df, int, on=:shots)
        end
        dropmissing!(df)
        coerce!(df, :label => OrderedFactor)
    end
end

# classification_v1!()
# classification_v2!()