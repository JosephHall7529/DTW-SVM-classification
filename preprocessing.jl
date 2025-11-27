# id codes 2D

function update_A2DF!(ts::Tuple{String, Int}, id::String)
        dict = OrderedDict("tok"=>ts[1], "shot"=>ts[2], "id"=>id)
        for feat in vcat(features_2D, "time_frames_".*["ITPA", "efit", "ecrh", "icrh", "nbi", "scal", "xcs", "icp", "tot", "fpg", "fpc", "ecs", "nis"])
                dict[feat] = 0
        end

        D = DATA_.data_2D.data[ts[1]].profiles[ts[2]]

        if typeof(D) == SAL_profile_extract
                diagnostics = D.diagnostics
        end

        if typeof(D) == AUG_profile_extract
                diagnostics = D.diagnostics
        end

        if typeof(D) == ITPA_profile_extract
                feats = D.features
                for feat in feats
                        println(ts[1], " ", ts[2], " ", feat)
                        dict[feat] = 1
                end
                feature_data = D.data[feats[1]]
                if length(feature_data.t) > 10
                        dict["time_frames_ITPA"] = 1
                end
        elseif typeof(D) == SAL_profile_extract
                diag_feats = D.features
                for diag in D.diagnostics
                        feats = diag_feats[diag]
                        for feat in feats
                                println("jet: ", ts[2], " $feat")
                                dict[feat] = 1
                        end
                        feature_data = D.data[diag][feats[1]]
                        if length(feature_data.t) > 10
                                dict["time_frames_"*diag] = 1
                        end
                end
        elseif typeof(D) == AUG_profile_extract
                diag_feats = D.features
                for diag in D.diagnostics
                        feats = diag_feats[diag]
                        for feat in feats
                                println("aug: ", ts[2], " $feat")
                                dict[feat] = 1
                        end
                        feature_data = D.data[diag][feats[1]]
                        if length(feature_data.t) > 10
                                dict["time_frames_"*diag] = 1
                        end
                end
        end
        append!(all_2D_features, DataFrame(dict))
end
function add_shot!(ts::Tuple{String, Int}, D::DataFrame=id_codes_2D)
        
        if ts[1] == "aug"
                search = ["AUG", "AUGW"]
        elseif ts[1] == "jet"
                search = ["JET", "JETILW"]
        else
                search = [uppercase(ts[1])]
        end

        # if !isempty(D)
        #         SS = (@subset D @byrow begin
        #                 in(:TOK, search)
        #                 :SHOT == ts[2]
        #                 end
        #         )

        #         @assert !isempty(SS) "shot already included"
        # end

        new_data = DataFrame(:TOK => ts[1], :SHOT => ts[2])
        new_data.id .= randstring(['A':'Z'; '0':'9'], 6)
        new_data.ITPA_data .= 0.0

        IDS = (@subset original_space(DATA_) @byrow begin
                in(:TOK, search)
                :SHOT == ts[2]
        end).id

        if !isempty(IDS)
                ID = IDS[1]
                new_data[1, :id] = ID
                new_data[1, :ITPA_data] = 1.0
        end

        append!(D, new_data)
        update_A2DF!(ts, new_data.id[1])
end

function files(str::AbstractString=pwd(), dotfiles::Bool=true; kwargs...)
        contents = readdir(str; kwargs...)
        if !dotfiles
                ind = findall(i -> i[1]!=='.', contents)
                return contents[ind]
        end
        return contents
end

# # creating id_codes_2D and all_2D_features
# all_2D_features = DataFrame()
# id_codes_2D = DataFrame()
# for tok in ["aug", "jet", "cmod", "d3d", "tftr", "txtr"]
#         shots = DATA_.data_2D.data[tok].shots
#         for shot in shots
#                 add_shot!((tok, shot))
#         end
# end
# id_codes_2D = select(unique(id_codes_2D, :id), [:id, :TOK, :SHOT], :)
# all_2D_features = select(unique(all_2D_features, :id), [:id, :tok, :shot], :)
# CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/id_codes_2D.csv", id_codes_2D)
# CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/common_2D_features.csv", all_2D_features)

# rename the dataframe columns, and remove the useless `index` column
# begin
#         for shot in AUG_shots
#                 for feat in ["IP", "BETAPOL", "Q95", "LI", "PNBI", "PECRH", "PICRH", "NGW"] 
#                         try 
#                                 D = CSV.read(AUG_profile_data*"/$(shot)/$(AUG_dict[feat])/$(feat).csv", DataFrame, stringtype=String)
#                                 println(shot, ": ", feat)
#                         catch 
#                                 continue
#                         end
#                         if feat == "PICRH"
#                                 try 
#                                         D = D[!, ["t", "PICRN"]]
#                                         rename!(D, ["t", "PICRH"])
#                                 catch 
#                                         D = D[!, ["t", "PICRH"]]
#                                 end
#                         elseif feat == "NGW"
#                                 try 
#                                         D = D[!, ["t", "n/nGW"]]
#                                         rename!(D, ["t", "NGW"]) 
#                                 catch
#                                         D = D[!, ["t", "NGW"]]
#                                 end
#                         end
#                         D = dropmissing(D[!, ["t", feat]])
#                         CSV.write(AUG_profile_data*"/$(shot)/$(AUG_dict[feat])/$(feat).csv", D)
#                 end
#         end
# end
# begin
#         for shot in JET_shots
#                 for feat in ["btnm", "btpm", "bttm", "f", "jsur", "li3m", "mi3m", "pio", "q95", "xip", "nbi_ptot", "icrh_ptot", "ecrh_ptot", "fgdl", "tril", "triu", "wp"] 
#                         try 
#                                 D = CSV.read(SAL_profile_data*"/$(shot)/$(SAL_dict[SAL_features_dict[feat]])/$(feat).csv", DataFrame, stringtype=String)
#                                 println(shot, ": ", feat)
#                         catch 
#                                 continue
#                         end
#                         if names(D) == ["t", SAL_features_dict[feat]]
#                                 continue
#                         end
#                         try 
#                                 D = D[!, Not("Column1")]
#                                 rename!(D, ["t", SAL_features_dict[feat]])
#                                 D = dropmissing(D[!, ["t", SAL_features_dict[feat]]])
#                                 CSV.write(SAL_profile_data*"/$(shot)/$(SAL_dict[SAL_features_dict[feat]])/$(feat).csv", D)
#                         catch
#                                 println("more than 2 columns")
#                         end
                        
#                 end
#         end
# end

# make sure every shot has a PECRH and PICRH csv file, even if its just a row of zeros.
# for shot in AUG_shots
#         icrh = isfile("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/icp/PICRH.csv")
#         ecrh = isfile("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/ecs/PECRH.csv")
#         if icrh && ecrh 
#                 continue
#         end
#         if icrh || ecrh 
#                 if icrh
#                         D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/fpg/Q95.csv", DataFrame)
#                         new = deepcopy(D)
#                         new[!, "Q95"] .= 0.0
#                         println(shot, ": ", size(new, 2))
#                         rename!(new, ["t", "PECRH"])
#                         CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/ecs/PECRH.csv", new)
#                 elseif ecrh
#                         D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/fpg/Q95.csv", DataFrame)
#                         new = deepcopy(D)
#                         new[!, "Q95"] .= 0.0
#                         rename!(new, ["t", "PICRH"])
#                         CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/icp/PICRH.csv", new)
#                 end
#         end
#         if !icrh && !ecrh 
#                 D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/fpg/Q95.csv", DataFrame)
#                 icrhD = deepcopy(D)
#                 icrhD.Q95 .= 0.0
#                 rename!(icrhD, ["t", "PICRH"])
#                 ecrhD = deepcopy(D)
#                 ecrhD.Q95 .= 0.0
#                 rename!(ecrhD, ["t", "PECRH"])
#                 CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/icp/PICRH.csv", icrhD)
#                 CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/ecs/PECRH.csv", ecrhD)
#         end
# end
# for shot in JET_shots
#         println(shot)
#         icrh = isfile("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/icrh/ptot.csv")
#         ecrh = isfile("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/ecrh/ptot.csv")
#         if icrh && ecrh 
#                 continue
#         end
#         if icrh || ecrh 
#                 if icrh
#                         if !isfile("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/efit/q95.csv")
#                                 continue
#                         end
#                         D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/efit/q95.csv", DataFrame)
#                         new = deepcopy(D)
#                         new[!, "Q95"] .= 0.0
#                         println(shot, ": ", size(new, 2))
#                         rename!(new, ["t", "PECRH"])
#                         CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/ecrh/ecrh_ptot.csv", new)
#                 elseif ecrh
#                         if !isfile("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/efit/q95.csv")
#                                 continue
#                         end
#                         D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/efit/q95.csv", DataFrame)
#                         new = deepcopy(D)
#                         new[!, "Q95"] .= 0.0
#                         rename!(new, ["t", "PICRH"])
#                         CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/icrh/icrh_ptot.csv", new)
#                 end
#         end
#         if !icrh && !ecrh 
#                 if !isfile("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/efit/q95.csv")
#                         continue
#                 end
#                 D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/efit/q95.csv", DataFrame)
#                 icrhD = deepcopy(D)
#                 icrhD.Q95 .= 0.0
#                 rename!(icrhD, ["t", "PICRH"])
#                 ecrhD = deepcopy(D)
#                 ecrhD.Q95 .= 0.0
#                 rename!(ecrhD, ["t", "PECRH"])
#                 CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/icrh/icrh_ptot.csv", icrhD)
#                 CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/JET/$(shot)/ecrh/ecrh_ptot.csv", ecrhD)
#         end
# end

## If the above code works but the csv files need more rows to be useful
# # for shot in AUG_shots
#         println(shot)
#         D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/ecs/PECRH.csv", DataFrame)
#         l = size(D, 1)
#         if zeros(l) == D.PECRH
#                 try 
#                         new = deepcopy(CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/fpg/Q95.csv", DataFrame))
#                 catch
#                         continue
#                 end
#                 new[!, "Q95"] .= 0.0
#                 rename!(new, ["t", "PECRH"])
#                 CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/ecs/PECRH.csv", new)
#         end
#         D = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/icp/PICRH.csv", DataFrame)
#         l = size(D, 1)
#         if zeros(l) == D.PICRH
#                 try
#                         new = deepcopy(CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/fpg/Q95.csv", DataFrame))
#                 catch 
#                         continue
#                 end
#                 new[!, "Q95"] .= 0.0
#                 rename!(new, ["t", "PICRH"])
#                 CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/AUG/$(shot)/icp/PICRH.csv", new)
#         end
# # end

struct unique_data
        corrected::DataFrame
        original::DataFrame
        identifiers::Vector
        duplicates::Vector
        function unique_data(df::DataFrame, identifiers::Vector; override=false, index_start=1)
                index!(df, override=override, start=index_start)
                reduced_by_identifiers = groupby(df, identifiers)
                identified = findall(i -> size(i, 1) > 1, reduced_by_identifiers)
                if isempty(identified)
                        return new(df, df, identifiers, Vector())
                end 
                duplicates = vcat([Pair(i.ind) for i in reduced_by_identifiers[identified]]...)
                rm_ind = [paired[1] for paired in duplicates]
                corrected = (@subset df @byrow !in(:ind, rm_ind))
                new(corrected, df, identifiers, duplicates)
        end
end
# Global database 
function unique_data!(df::DataFrame, identifiers::Vector; override=false, index_start=1)
        index!(df, override=override, start=index_start)
        indexstring!(df, override=override)
        reduced_by_identifiers = groupby(df, identifiers)
        identified = findall(i -> size(i, 1) > 1, reduced_by_identifiers)
        duplicates = vcat([i.ind |> Pair for i in reduced_by_identifiers[identified]]...)
        rm_ind = [paired[2] for paired in duplicates]
        return @subset! df @byrow !in(:ind, rm_ind)
end

import Base.log
function log(df::DataFrame, features::Vector=log_features)
    for feat in features
        df[!, feat] = log.(abs.(df[!, feat]))
    end
    return df[!, vcat(["id", "TOK", "SHOT"], String.(features), ["HYBRID"])]
end