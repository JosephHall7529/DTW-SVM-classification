# iterate over all trials: perform hyper-parameter search for set variables (add into the code manually)
begin
    initialize_trials = true
    if initialize_trials
        complete_CH_trial_1 = Dict()
        complete_CH_trial_2 = Dict()
        complete_CH_trial_3 = Dict()
        complete_CH_trial_4 = Dict()
        complete_CH_trial_5 = Dict()
    end
    for trial in 1:1
        begin
            if trial == 1
                CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots
            
                COchoices = sort(tr_1[1:COn])
                CEchoices = sort(tr_1[COn+1:COn+CEn])
                CLchoices = sort(tr_1[COn+CEn+1:end])
            elseif trial == 2
                CLn, CEn, COn = trial_dataset_meta[:trial_2][2:4]
                tr_2 = (@subset trials_D @byrow :trial_2 == 1).shots
            
                COchoices = sort(tr_2[1:COn])
                CEchoices = sort(tr_2[COn+1:COn+CEn])
                CLchoices = sort(tr_2[COn+CEn+1:end])
            elseif trial == 3
                CLn, CEn, COn = trial_dataset_meta[:trial_3][2:4]
                tr_3 = (@subset trials_D @byrow :trial_3 == 1).shots
            
                COchoices = sort(tr_3[1:COn])
                CEchoices = sort(tr_3[COn+1:COn+CEn])
                CLchoices = sort(tr_3[COn+CEn+1:end])
            elseif trial == 4
                CLn, CEn, COn = trial_dataset_meta[:trial_4][2:4]
                tr_4 = (@subset trials_D @byrow :trial_4 == 1).shots
            
                COchoices = sort(tr_4[1:COn])
                CEchoices = sort(tr_4[COn+1:COn+CEn])
                CLchoices = sort(tr_4[COn+CEn+1:end])
            elseif trial == 5
                CLn, CEn, COn = trial_dataset_meta[:trial_5][2:4]
                tr_5 = (@subset trials_D @byrow :trial_5 == 1).shots
            
                COchoices = sort(tr_5[1:COn])
                CEchoices = sort(tr_5[COn+1:COn+CEn])
                CLchoices = sort(tr_5[COn+CEn+1:end])
            end
        end
        labelled_data = let
            labelled_shots = vcat(CLchoices, CEchoices, COchoices)
            labels = vcat(repeat(["C-LH"], length(CLchoices)),
                repeat(["C-EH"], length(CEchoices)),
                repeat(["CO"], length(COchoices))
            ) 
            d = OrderedDict([ts for ts in labelled_shots] .=> labels)
        end 
        comp_CH = emp("complete_CH_trial_$trial")
        
        # ADD PARAMETERS HERE
        # params = collect(powerset(["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "LI", "NGW"], 1, 1))
        # begin
        #     if in(trial, [3])
        #         params = [vcat(["NGW"], [para]) for para in ["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "LI"]]
        #     elseif in(trial, [1, 2, 4, 5])
        #         params = [vcat(["BETAPOL"], [para]) for para in ["IP", "PNBI", "PECRH", "PICRH", "Q95", "LI", "NGW"]]
        #     end
        # end
        # begin
        #     params = [vcat(["PNBI", "PICRH"], [para]) for para in ["IP", "PECRH", "Q95", "LI", "NGW"]]
        # end
        params = [["PNBI"]]
        for para in params
            para = naming(para)
            if haskey(data_dtw, para)
                data_int = data_dtw[para]
            else
                println("data_dtw doesnt have ", para)
                continue
            end
            if !isempty(comp_CH)
                if in(para, collect(keys(comp_CH)))
                    continue
                end
            end
            println("\n", para)
            k = [Int(div(CLn, 10/7)), Int(div(CEn, 10/7)), Int(div(COn, 10/7))]
            println(k)
            levels = ["C-LH", "C-EH", "CO"]
            
            comp_CH[para] = try 
                hyper_parameter_search(data_int, labelled_data, k, interesting="CO", N=60, metric="MulticlassFalsePositiveRate(;levels=$levels, perm=[1, 2, 3], checks=false)", max=false)
            catch
            end
        end
    end

    CH = complete_CH_trial_1[naming(["PNBI"])]
    data_int = data_dtw[naming(["PNBI"])]

    tss = tok_shots((@subset which(["PNBI"]) @byrow in(:tok, ["aug"])))
    ℓT = length(tss)
    dict_shot = Dict(1:ℓT .=> tss)

    begin
        CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                    tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots
                
                    COchoices = sort(tr_1[1:COn])
                    CEchoices = sort(tr_1[COn+1:COn+CEn])
                    CLchoices = sort(tr_1[COn+CEn+1:end])
        labelled_data = let
            labelled_shots = vcat(CLchoices, CEchoices, COchoices)
            labels = vcat(repeat(["C-LH"], length(CLchoices)),
                repeat(["C-EH"], length(CEchoices)),
                repeat(["CO"], length(COchoices))
            ) 
            d = OrderedDict([ts for ts in labelled_shots] .=> labels)
        end 
    end
    begin
        labelled_ts = collect(labelled_data |> keys)
        labelled_y = collect(labelled_data |> values)
        shot_dict = Dict([a => n for (n, a) in enumerate(data_int.tok_shots)])
        labelled_ind = [shot_dict[k] for k in labelled_ts]
    end
end

import MLJParticleSwarmOptimization
struct classification_performance
    features::Vector{String}
    trial::Symbol
    best_model::DTW_SVM
    train_results::DataFrame
    test_results::DataFrame
end

trial_one = Dict{Symbol, Dict{String, classification_performance}}()
# for tr_name in [:trial_1, :trial_2, :trial_3, :trial_4, :trial_5]
    # result_tmp = Dict()
    # for para in collect(powerset(["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "LI", "NGW"], 1, 1))
        # if in(tr_name, keys(trial_one))
        #     if in(naming(para), keys(trial_one[tr_name]))
        #         println(para)
        #         continue
        #     end
        # end
        tr_name=:trial_2
        para=["IP"]
        labelled = (@subset trials_D @byrow $tr_name == 1)
        unseen = (@subset trials_D @byrow $tr_name !== 1)

        if in("NGW", para)
            rm_ind = findfirst(i -> i == ("aug", 29624), unseen.shots)
            deleteat!(unseen, rm_ind)
        end

        model_new = half_model(para, labelled.shots, labelled.label;
            target_shot_length=100,
            transportcost=1.1)
        begin
            c0 = quantile(Array(model_new.database.training_data.cosine_cost[:, 2:end])[:], 0.05)
            c1 = quantile(Array(model_new.database.training_data.cosine_cost[:, 2:end])[:], 0.95)
            ft0 = quantile(Array(model_new.database.training_data.flat_top_cost[:, 2:end])[:], 0.05)
            ft1 = quantile(Array(model_new.database.training_data.flat_top_cost[:, 2:end])[:], 0.95)

            # r1 = range(model_tmp, :transportcost, lower=1, upper=1.5)
            # r2 = range(model_tmp, :target_shot_length, lower=100, upper=1000)
            r2 = range(model_new, :C_cosine, lower=c0, upper=c1)
            r3 = range(model_new, :C_flat_top, lower=ft0, upper=ft1)
            r4 = range(model_new, :C_cost, lower=0.5, upper=5) 

            # n_particles = 100_000 is a lot
            self_tuning_tree = TunedModel(
                    model=model_new,
                    tuning=MLJParticleSwarmOptimization.AdaptiveParticleSwarm(n_particles=30),
                    # tuning=Grid(resolution=10),
                    resampling=StratifiedCV(nfolds=5, shuffle=true),
                    range=[r2, r2, r4],
                    measure=MulticlassFalsePositiveRate(levels=levels),
                    n=10,
                    check_measure=true
            )

            mach = machine(self_tuning_tree, labelled.shots, labelled.label) 
            # try 
                MLJBase.fit!(mach, verbosity=1)
            # catch
            #     continue
            # end
        end
        mach_save = deepcopy(mach)
        hyp_new = deepcopy(fitted_params(mach).best_model)
        
        begin
            fitresult, cache, report = MLJBase.fit(fitted_params(mach).best_model, 0, labelled.shots, labelled.label)
            ỹ = MLJBase.predict(fitted_params(mach).best_model, fitresult, unseen.shots)

            mcFPR, BA, A = MulticlassFalsePositiveRate(levels=levels)(ỹ, unseen.label), BalancedAccuracy()(ỹ, unseen.label), Accuracy()(ỹ, unseen.label)
        end

        test_result = DataFrame([mcFPR BA A], [:mcFPR, :BA, :ACC])
        result_tmp[naming(para)] = classification_performance(para, tr_name, fitted_params(mach).best_model, train_result, test_result) 
#     end
#     trial_one[tr_name] = result_tmp
# end
trial_one[:trial_2]["IP"]
begin
    f = Figure();
    a = Axis(f[1, 1])
    for (n, (feat, D)) in enumerate(trial_1_results)
        hist!(a, D.train_results[:, 4], label=feat, normalization=:probability, color=(colors_15[2n], 0.55))
    end
    f[1, 2] = Legend(f, a, "features", framevisible=false)
    f
end

trial_2_results = Dict()
for para in vcat.(["BETAPOL"], collect(powerset(["IP", "PNBI", "PECRH", "PICRH", "Q95", "LI", "NGW"], 1, 1)))
    if in(naming(para), keys(trial_2_results))
        println(para)
        continue
    end
    labelled = (@subset trials_D @byrow :trial_1 == 1)
    unseen = (@subset trials_D @byrow :trial_1 !== 1)

    if in("NGW", para)
        rm_ind = findfirst(i -> i == ("aug", 29624), unseen.shots)
        deleteat!(unseen, rm_ind)
    end

    model_new = half_model(para, labelled.shots, labelled.label;
    target_shot_length=100,
    transportcost=1.1)
    begin
        c0 = quantile(Array(model_new.database.training_data.cosine_cost[:, 2:end])[:], 0.05)
        c1 = quantile(Array(model_new.database.training_data.cosine_cost[:, 2:end])[:], 0.95)
        ft0 = quantile(Array(model_new.database.training_data.flat_top_cost[:, 2:end])[:], 0.05)
        ft1 = quantile(Array(model_new.database.training_data.flat_top_cost[:, 2:end])[:], 0.95)

        # r1 = range(model_tmp, :transportcost, lower=1, upper=1.5)
        # r2 = range(model_tmp, :target_shot_length, lower=100, upper=1000)
        r2 = range(model_new, :C_cosine, lower=c0, upper=c1)
        r3 = range(model_new, :C_flat_top, lower=ft0, upper=ft1)
        r4 = range(model_new, :C_cost, lower=0.5, upper=5) 

        # n_particles = 100_000 is a lot
        self_tuning_tree = TunedModel(
                model=model_new,
                tuning=MLJParticleSwarmOptimization.AdaptiveParticleSwarm(n_particles=30),
                # tuning=Grid(resolution=10),
                resampling=StratifiedCV(nfolds=5, shuffle=true),
                range=[r2, r2, r4],
                measure=MulticlassFalsePositiveRate(levels=levels),
                n=10,
                check_measure=true
        )

        mach = machine(self_tuning_tree, labelled.shots, labelled.label) 
        try 
            MLJBase.fit!(mach, verbosity=1)
        catch
            continue
        end
    end
    mach_save = deepcopy(mach)
    hyp_new = deepcopy(fitted_params(mach).best_model)

    train_result = evaluate(hyp_new, [MulticlassFalsePositiveRate(levels=levels), BalancedAccuracy(), Accuracy()])
    
    begin
        fitresult, cache, report = MLJBase.fit(fitted_params(mach).best_model, 0, labelled.shots, labelled.label)
        ỹ = MLJBase.predict(fitted_params(mach).best_model, fitresult, unseen.shots)

        mcFPR, BA, A = MulticlassFalsePositiveRate(levels=levels)(ỹ, unseen.label), BalancedAccuracy()(ỹ, unseen.label), Accuracy()(ỹ, unseen.label)
    end

    test_result = DataFrame([mcFPR BA A], [:mcFPR, :BA, :ACC])
    trial_1_results[naming(para)] = classification_performance(para, :trial_1, fitted_params(mach).best_model, train_result, test_result) 
end



let 
    # f = Figure();
    # a = Axis(f[1, 1])

    tr_1, te_1 = partition(labelled, 0.7; stratify=labelled.label)
    model = fitted_params(mach).best_model
    svmmodel = mach.fitresult.fitresult
    # # overview_fit(model, svmmodel, unseen_1.shots, unseen_1.label)
    overview_fit(model, svmmodel, tr_1.shots, tr_1.label)
    # signal_1 = 
    # mat = dtw_cost_matrix(signal_2, signal_1, Cityblock())
	# cost, i1, i2 = DynamicAxisWarping.trackback(mat)
end


fitted_params(mach).best_model
# using known hyp-fit-values
begin
    model_check = half_model(["PNBI"], labelled.shots, labelled.label;
                target_shot_length = 160, 
                transportcost = 1.267144391963193, 
                C_cosine = 0.08891198597494063, 
                C_flat_top = 5.742405215528001, 
                C_cost = 2.1018886822814964, 
                short_IP = false,
                kernel=false
    )
    Array(model_check.database.training_data.cosine_cost[!, 2:end]) |> issymmetric

    # tr_1, te_1 = partition(labelled, 0.7;
    #                         shuffle=true,
    #                         stratify=labelled.label)

    # fitresult, cache, report = MLJBase.fit(model_check, 0, labelled.shots, labelled.label)
    # ỹ = MLJBase.predict(model_check, fitresult, unseen_1.shots)

    # MulticlassFalsePositiveRate(levels=["LH", "EH", "CO"])(ỹ, unseen_1.label), BalancedAccuracy()(ỹ, unseen_1.label)
end



# Table of test trial results
begin
    Results = DataFrame()

    run_trial_single_param = false 
    if run_trial_single_param

        Results = DataFrame()
        params = collect(powerset(["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "LI", "NGW"], 1, 1))
        for ((n, tr_name), predictor) in Iterators.product(zip(1:5, [:trial_1, :trial_2, :trial_3, :trial_4, :trial_5]), params)
        # tr_name = :trial_1
        # n = 1
        # predictor = ["BETAPOL"]
            para = naming(predictor)
            comp_CH = emp("complete_CH_$(String(tr_name))_")
            
            if !in(para, keys(comp_CH))
                # continue
            end
            
            trial = n
            begin
                if trial == 1
                    CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                    tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots

                    COchoices = sort(tr_1[1:COn])
                    CEchoices = sort(tr_1[COn+1:COn+CEn])
                    CLchoices = sort(tr_1[COn+CEn+1:end])
                elseif trial == 2
                    CLn, CEn, COn = trial_dataset_meta[:trial_2][2:4]
                    tr_2 = (@subset trials_D @byrow :trial_2 == 1).shots

                    COchoices = sort(tr_2[1:COn])
                    CEchoices = sort(tr_2[COn+1:COn+CEn])
                    CLchoices = sort(tr_2[COn+CEn+1:end])
                elseif trial == 3
                    CLn, CEn, COn = trial_dataset_meta[:trial_3][2:4]
                    tr_3 = (@subset trials_D @byrow :trial_3 == 1).shots

                    COchoices = sort(tr_3[1:COn])
                    CEchoices = sort(tr_3[COn+1:COn+CEn])
                    CLchoices = sort(tr_3[COn+CEn+1:end])
                elseif trial == 4
                    CLn, CEn, COn = trial_dataset_meta[:trial_4][2:4]
                    tr_4 = (@subset trials_D @byrow :trial_4 == 1).shots

                    COchoices = sort(tr_4[1:COn])
                    CEchoices = sort(tr_4[COn+1:COn+CEn])
                    CLchoices = sort(tr_4[COn+CEn+1:end])
                elseif trial == 5
                    CLn, CEn, COn = trial_dataset_meta[:trial_5][2:4]
                    tr_5 = (@subset trials_D @byrow :trial_5 == 1).shots

                    COchoices = sort(tr_5[1:COn])
                    CEchoices = sort(tr_5[COn+1:COn+CEn])
                    CLchoices = sort(tr_5[COn+CEn+1:end])
                end
                labelled_data = let
                    labelled_shots = vcat(CLchoices, CEchoices, COchoices)
                    labels = vcat(repeat(["C-LH"], length(CLchoices)),
                        repeat(["C-EH"], length(CEchoices)),
                        repeat(["CO"], length(COchoices))
                    ) 
                    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
                end
            end
            
            begin
                if haskey(comp_CH, para)
                    # if isempty(comp_CH[para].hyperparameters_meta)
                    if comp_CH[para] === missing
                        continue
                    end
                     hyperparameters_meta!(comp_CH[para]; N=1_000)
                    # end
                    dscrb = describe(comp_CH[para].hyperparameters_meta)
                    dict_min = Dict(dscrb.variable .=> dscrb.min)
                    dict_med = Dict(dscrb.variable .=> dscrb.median)
                    dict_mean = Dict(dscrb.variable .=> round.(dscrb.mean, digits=3))
                    dict_max =  Dict(dscrb.variable .=> dscrb.max)
                    nSV_EH = dict_med[Symbol("#SV EH")]
                    nSV_LH = dict_med[Symbol("#SV LH")]
                    nSV_CO = dict_med[Symbol("#SV CO")]

                    n_EH = dict_med[Symbol("#EH")]
                    n_LH = dict_med[Symbol("#LH")]
                    n_CO = dict_med[Symbol("#CO")]

                    acc = (dict_med[Symbol("Accuracy")], dict_min[Symbol("Accuracy")], dict_med[Symbol("Accuracy")], dict_max[Symbol("Accuracy")])
                    b_acc =  (dict_med[Symbol("Balanced Accuracy")], dict_min[Symbol("Balanced Accuracy")], dict_med[Symbol("Balanced Accuracy")], dict_max[Symbol("Balanced Accuracy")])
                    fpr = (dict_med[Symbol("False Positive Rate")], dict_min[Symbol("False Positive Rate")], dict_med[Symbol("False Positive Rate")], dict_max[Symbol("False Positive Rate")])
                    # b_acc_CO_EH = (dict_med[Symbol("Balanced Accuracy CO EH")], dict_min[Symbol("Balanced Accuracy CO EH")], dict_med[Symbol("Balanced Accuracy CO EH")], dict_max[Symbol("Balanced Accuracy CO EH")])

                    IP = in("IP", predictor) ? 1 : 0
                    PNBI = in("PNBI", predictor) ? 1 : 0
                    NGW = in("NGW", predictor) ? 1 : 0
                    PECRH = in("PECRH", predictor) ? 1 : 0
                    PICRH = in("PICRH", predictor) ? 1 : 0
                    Q95 = in("Q95", predictor) ? 1 : 0
                    LI = in("LI", predictor) ? 1 : 0
                    BETAPOL = in("BETAPOL", predictor) ? 1 : 0

                    int = OrderedDict("trial" => tr_name, L"I_p" => IP, L"P_\mathrm{NBI}" => PNBI, L"P_\mathrm{ECRH}" => PECRH, L"P_\mathrm{ICRH}" => PICRH, L"q_{95}" => Q95, L"\beta_\mathrm{p}" => BETAPOL, L"f_\mathrm{GW}" => NGW, L"\ell_i" => LI, 
                        "Accuracy" => acc, "Balanced Accuracy" => b_acc, "False Positive Rate" => fpr,
                        "#SV LH" => nSV_LH, "#SV EH" => nSV_EH, "#SV CO" => nSV_CO, 
                        "#LH" => n_LH, "#EH" => n_EH, "#CO" => n_CO)
                    append!(Results, DataFrame(int))
                end
            end
        end
    end
    # CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/metadata/trial_case/feature_selection/single_feature_run.csv", Results)
    
    run_trial_two_param = true 
    if run_trial_two_param
        Results = DataFrame()
        # params_POL = [vcat(["BETAPOL"], [para]) for para in ["IP", "PNBI", "PECRH", "PICRH", "Q95", "LI", "NGW"]]
        # params_NGW = [vcat(["NGW"], [para]) for para in ["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "LI"]]
        params = [["IP", "PNBI"]]
        for (n, (tr_name, predictors_)) in enumerate(zip([:trial_1, :trial_2, :trial_3, :trial_4, :trial_5], [params, params, params, params, params]))
            # println(tr_name, ": $n  ", predictors_)
            for predictors in predictors_ 
                para, S = naming(predictors), 123
                comp_CH = emp("complete_CH_$(String(tr_name))")
                
                # if !in(para, keys(comp_CH))
                #     continue
                # end
                
                trial = n
                begin
                    if trial == 1
                        CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                        tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots

                        COchoices = sort(tr_1[1:COn])
                        CEchoices = sort(tr_1[COn+1:COn+CEn])
                        CLchoices = sort(tr_1[COn+CEn+1:end])
                    elseif trial == 2
                        CLn, CEn, COn = trial_dataset_meta[:trial_2][2:4]
                        tr_2 = (@subset trials_D @byrow :trial_2 == 1).shots

                        COchoices = sort(tr_2[1:COn])
                        CEchoices = sort(tr_2[COn+1:COn+CEn])
                        CLchoices = sort(tr_2[COn+CEn+1:end])
                    elseif trial == 3
                        CLn, CEn, COn = trial_dataset_meta[:trial_3][2:4]
                        tr_3 = (@subset trials_D @byrow :trial_3 == 1).shots

                        COchoices = sort(tr_3[1:COn])
                        CEchoices = sort(tr_3[COn+1:COn+CEn])
                        CLchoices = sort(tr_3[COn+CEn+1:end])
                    elseif trial == 4
                        CLn, CEn, COn = trial_dataset_meta[:trial_4][2:4]
                        tr_4 = (@subset trials_D @byrow :trial_4 == 1).shots

                        COchoices = sort(tr_4[1:COn])
                        CEchoices = sort(tr_4[COn+1:COn+CEn])
                        CLchoices = sort(tr_4[COn+CEn+1:end])
                    elseif trial == 5
                        CLn, CEn, COn = trial_dataset_meta[:trial_5][2:4]
                        tr_5 = (@subset trials_D @byrow :trial_5 == 1).shots

                        COchoices = sort(tr_5[1:COn])
                        CEchoices = sort(tr_5[COn+1:COn+CEn])
                        CLchoices = sort(tr_5[COn+CEn+1:end])
                    end
                    labelled_data = let
                        labelled_shots = vcat(CLchoices, CEchoices, COchoices)
                        labels = vcat(repeat(["C-LH"], length(CLchoices)),
                            repeat(["C-EH"], length(CEchoices)),
                            repeat(["CO"], length(COchoices))
                        ) 
                        d = OrderedDict([ts for ts in labelled_shots] .=> labels)
                    end
                end
                
                begin
                    # if isempty(comp_CH[para].hyperparameters_meta)
                        hyperparameters_meta!(comp_CH[para]; N=1000)
                    # end
                    dscrb = describe(comp_CH[para].hyperparameters_meta)
                    dict_min = Dict(dscrb.variable .=> dscrb.min)
                    dict_med = Dict(dscrb.variable .=> dscrb.median)
                    dict_mean = Dict(dscrb.variable .=> round.(dscrb.mean, digits=3))
                    dict_max =  Dict(dscrb.variable .=> dscrb.max)
                    nSV_EH = dict_med[Symbol("#SV EH")]
                    nSV_LH = dict_med[Symbol("#SV LH")]
                    nSV_CO = dict_med[Symbol("#SV CO")]

                    n_EH = dict_med[Symbol("#EH")]
                    n_LH = dict_med[Symbol("#LH")]
                    n_CO = dict_med[Symbol("#CO")]

                    acc = (dict_med[Symbol("Accuracy")], dict_min[Symbol("Accuracy")], dict_med[Symbol("Accuracy")], dict_max[Symbol("Accuracy")])
                    b_acc =  (dict_med[Symbol("Balanced Accuracy")], dict_min[Symbol("Balanced Accuracy")], dict_med[Symbol("Balanced Accuracy")], dict_max[Symbol("Balanced Accuracy")])
                    fpr = (dict_med[Symbol("False Positive Rate")], dict_min[Symbol("False Positive Rate")], dict_med[Symbol("False Positive Rate")], dict_max[Symbol("False Positive Rate")])
                    # b_acc_CO_EH = (dict_med[Symbol("Balanced Accuracy CO EH")], dict_min[Symbol("Balanced Accuracy CO EH")], dict_med[Symbol("Balanced Accuracy CO EH")], dict_max[Symbol("Balanced Accuracy CO EH")])

                    IP = in("IP", predictors) ? 1 : 0
                    PNBI = in("PNBI", predictors) ? 1 : 0
                    NGW = in("NGW", predictors) ? 1 : 0
                    PECRH = in("PECRH", predictors) ? 1 : 0
                    PICRH = in("PICRH", predictors) ? 1 : 0
                    Q95 = in("Q95", predictors) ? 1 : 0
                    LI = in("LI", predictors) ? 1 : 0
                    BETAPOL = in("BETAPOL", predictors) ? 1 : 0

                    int = OrderedDict("trial" => tr_name, L"I_p" => IP, L"P_\mathrm{NBI}" => PNBI, L"P_\mathrm{ECRH}" => PECRH, L"P_\mathrm{ICRH}" => PICRH, L"q_{95}" => Q95, L"\beta_\mathrm{p}" => BETAPOL, L"f_\mathrm{GW}" => NGW, L"\ell_i" => LI, 
                        "Accuracy" => acc, "Balanced Accuracy" => b_acc, "False Positive Rate" => fpr,
                        "#SV LH" => nSV_LH, "#SV EH" => nSV_EH, "#SV CO" => nSV_CO, 
                        "#LH" => n_LH, "#EH" => n_EH, "#CO" => n_CO)
                    append!(Results, DataFrame(int))
                end
            end
        end
        Results
    end
    # CSV.wrinte("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/metadata/trial_case/feature_selection/two_feature_run.csv", Results)

    run_trial_three_param = false 
    if run_trial_three_param
        Results = DataFrame()
        # params_IP = [vcat(["IP"], [para]) for para in ["PNBI", "PECRH", "PICRH", "Q95", "LI", "NGW"]]
        params = [vcat(["PNBI", "PICRH"], [para]) for para in ["IP", "PECRH", "Q95", "LI", "NGW"]]
        for ((n, tr_name, predictors_)) in (zip([2, 5], [:trial_2, :trial_5], [params, params]))
            # println(tr_name, ": $n  ", predictors_)
            for predictors in predictors_ 
                para = naming(predictors)
                comp_CH = emp("complete_CH_$(String(tr_name))")
                
                if !in(para, keys(comp_CH))
                    # continue
                end
                
                trial = n
                begin
                    if trial == 1
                        CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                        tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots

                        COchoices = sort(tr_1[1:COn])
                        CEchoices = sort(tr_1[COn+1:COn+CEn])
                        CLchoices = sort(tr_1[COn+CEn+1:end])
                    elseif trial == 2
                        CLn, CEn, COn = trial_dataset_meta[:trial_2][2:4]
                        tr_2 = (@subset trials_D @byrow :trial_2 == 1).shots

                        COchoices = sort(tr_2[1:COn])
                        CEchoices = sort(tr_2[COn+1:COn+CEn])
                        CLchoices = sort(tr_2[COn+CEn+1:end])
                    elseif trial == 3
                        CLn, CEn, COn = trial_dataset_meta[:trial_3][2:4]
                        tr_3 = (@subset trials_D @byrow :trial_3 == 1).shots

                        COchoices = sort(tr_3[1:COn])
                        CEchoices = sort(tr_3[COn+1:COn+CEn])
                        CLchoices = sort(tr_3[COn+CEn+1:end])
                    elseif trial == 4
                        CLn, CEn, COn = trial_dataset_meta[:trial_4][2:4]
                        tr_4 = (@subset trials_D @byrow :trial_4 == 1).shots

                        COchoices = sort(tr_4[1:COn])
                        CEchoices = sort(tr_4[COn+1:COn+CEn])
                        CLchoices = sort(tr_4[COn+CEn+1:end])
                    elseif trial == 5
                        CLn, CEn, COn = trial_dataset_meta[:trial_5][2:4]
                        tr_5 = (@subset trials_D @byrow :trial_5 == 1).shots

                        COchoices = sort(tr_5[1:COn])
                        CEchoices = sort(tr_5[COn+1:COn+CEn])
                        CLchoices = sort(tr_5[COn+CEn+1:end])
                    end
                    labelled_data = let
                        labelled_shots = vcat(CLchoices, CEchoices, COchoices)
                        labels = vcat(repeat(["C-LH"], length(CLchoices)),
                            repeat(["C-EH"], length(CEchoices)),
                            repeat(["CO"], length(COchoices))
                        ) 
                        d = OrderedDict([ts for ts in labelled_shots] .=> labels)
                    end
                end
                
                begin
                    # if isempty(comp_CH[para].hyperparameters_meta)
                        hyperparameters_meta!(comp_CH[para], labelled_data; N=1000)
                    # end
                    dscrb = describe(comp_CH[para].hyperparameters_meta) 
                    dict_med = Dict(dscrb.variable .=> dscrb.:median)
                    dict_mean = Dict(dscrb.variable .=> dscrb.:mean)
                    nSV_EH = dict_med[Symbol("#SV EH")]
                    nSV_LH = dict_med[Symbol("#SV LH")]
                    nSV_CO = dict_med[Symbol("#SV CO")]

                    n_EH = dict_med[Symbol("#EH")]
                    n_LH = dict_med[Symbol("#LH")]
                    n_CO = dict_med[Symbol("#CO")]

                    acc = dict_mean[Symbol("Accuracy")]
                    b_acc = dict_mean[Symbol("Balanced Accuracy")]
                    fpr = dict_mean[Symbol("False Positive Rate")]

                    IP = in("IP", predictors) ? 1 : 0
                    PNBI = in("PNBI", predictors) ? 1 : 0
                    NGW = in("NGW", predictors) ? 1 : 0
                    PECRH = in("PECRH", predictors) ? 1 : 0
                    PICRH = in("PICRH", predictors) ? 1 : 0
                    Q95 = in("Q95", predictors) ? 1 : 0
                    LI = in("LI", predictors) ? 1 : 0
                    BETAPOL = in("BETAPOL", predictors) ? 1 : 0

                    int = OrderedDict("trial" => tr_name, L"I_p" => IP, L"P_\mathrm{NBI}" => PNBI, L"P_\mathrm{ECRH}" => PECRH, L"P_\mathrm{ICRH}" => PICRH, L"q_{95}" => Q95, L"\beta_\mathrm{p}" => BETAPOL, L"f_\mathrm{GW}" => NGW, L"\ell_i" => LI, 
                        "Accuracy" => acc, "Balanced Accuracy" => b_acc, "False Positive Rate" => fpr,
                        "#SV LH" => nSV_LH, "#SV EH" => nSV_EH, "#SV CO" => nSV_CO, 
                        "#LH" => n_LH, "#EH" => n_EH, "#CO" => n_CO)
                    append!(Results, DataFrame(int))
                end
            end
        end
        Results
    end
    # CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/metadata/trial_case/feature_selection/three_feature_run.csv", Results) 

    load_trial = false    
    if load_trial
        Results = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/metadata/trial_case/feature_selection/single_feature_run.csv", DataFrame, stringtype=String)
        # Results = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/metadata/trial_case/feature_selection/two_feature_run.csv", DataFrame, stringtype=String)
        # Results = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/metadata/trial_case/feature_selection/three_feature_run.csv", DataFrame, stringtype=String)
    end

    sort!(Results, ["trial", "Balanced Accuracy"])
    Results |> vscodedisplay
    # # Results[[5, 12, 17, 24, 30], :] |> clipboard
end 

# Plot the train test case
let 
	# CairoMakie.activate!()
	f = Figure(size=(1000, 900))
	tests = [15714, 33371]
	opaque = ["IP", "PNBI"]
	σ = false
	begin 
        trial = 2
        para = naming(["IP", "PNBI"])
        S=1
        begin
            if trial == 1
                CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots

                COchoices = sort(tr_1[1:COn])
                CEchoices = sort(tr_1[COn+1:COn+CEn])
                CLchoices = sort(tr_1[COn+CEn+1:end])
            elseif trial == 2
                CLn, CEn, COn = trial_dataset_meta[:trial_2][2:4]
                tr_2 = (@subset trials_D @byrow :trial_2 == 1).shots

                COchoices = sort(tr_2[1:COn])
                CEchoices = sort(tr_2[COn+1:COn+CEn])
                CLchoices = sort(tr_2[COn+CEn+1:end])
            elseif trial == 3
                CLn, CEn, COn = trial_dataset_meta[:trial_3][2:4]
                tr_3 = (@subset trials_D @byrow :trial_3 == 1).shots

                COchoices = sort(tr_3[1:COn])
                CEchoices = sort(tr_3[COn+1:COn+CEn])
                CLchoices = sort(tr_3[COn+CEn+1:end])
            elseif trial == 4
                CLn, CEn, COn = trial_dataset_meta[:trial_4][2:4]
                tr_4 = (@subset trials_D @byrow :trial_4 == 1).shots

                COchoices = sort(tr_4[1:COn])
                CEchoices = sort(tr_4[COn+1:COn+CEn])
                CLchoices = sort(tr_4[COn+CEn+1:end])
            elseif trial == 5
                CLn, CEn, COn = trial_dataset_meta[:trial_5][2:4]
                tr_5 = (@subset trials_D @byrow :trial_5 == 1).shots

                COchoices = sort(tr_5[1:COn])
                CEchoices = sort(tr_5[COn+1:COn+CEn])
                CLchoices = sort(tr_5[COn+CEn+1:end])
            end
            labelled_data = let
                labelled_shots = vcat(CLchoices, CEchoices, COchoices)
                labels = vcat(repeat(["C-LH"], length(CLchoices)),
                    repeat(["C-EH"], length(CEchoices)),
                    repeat(["CO"], length(COchoices))
                ) 
                d = OrderedDict([ts for ts in labelled_shots] .=> labels)
            end
        end
        comp_CH = emp("complete_CH_trial_$(trial)")
		k = [Int(div(*(CLn, 7), 10)), Int(div(*(CEn, 7), 10)), Int(div(*(COn, 7),10))]
		train_tss = vcat(CLchoices, CEchoices, COchoices)
		n = length(train_tss)
    end
    begin
        cos_med, ft_med, C = comp_CH[para].hyperparameters
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		shot_dict = Dict([a => n for (n, a) in enumerate(comp_CH[para].data.tok_shots)])
		labelled_ind = [shot_dict[k] for k in labelled_ts]
		train_ind = training_partion(labelled_data, unique(labelled_y), k, S=S)

		X = exp.(-(Array(comp_CH[para].data.cosine_cost[!, 2:end]) ./ cos_med).^2) .*
			exp.(-(Array(comp_CH[para].data.flat_top_cost[!, 2:end]) ./ ft_med).^2)

		K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
		model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C)

		col_dict = Dict("CO"=>1, "C-EH"=>2, "C-LH"=>3)
		feat_norm_dict = Dict("IP" => 1, "PNBI" => 10, "PECRH" => 10, "PICRH" => 10, 
			"Q95" => 10, "NGW" => 1, "LI" => 1, "BETAPOL" => 1
		)
		feat_color_dict = Dict(
			"IP" => (colors[1], in("IP", opaque) ? 1 : 0.3),
			"PNBI" => (colors[2], in("PNBI", opaque) ? 1 : 0.3),
			"PECRH" => (colors[3], in("PECRH", opaque) ? 1 : 0.3),
			"PICRH" => (colors[4], in("PICRH", opaque) ? 1 : 0.3),
			"Q95" => (colors[5], in("Q95", opaque) ? 1 : 0.3),
			"NGW" => (colors[6], in("NGW", opaque) ? 1 : 0.3),
			"LI" => (colors[7], in("LI", opaque) ? 1 : 0.3),
			"BETAPOL" => (colors[8], in("BETAPOL", opaque) ? 1 : 0.3)
		)
		feat_style_dict = Dict(
			"IP" => :solid, 
			"PNBI" => :solid, 
			"PECRH" => :solid, 
			"PICRH" => :solid, 
			"Q95" => :dash, 
			"NGW" => :dash, 
			"LI" => :dash, 
			"BETAPOL" => :dash
		)

		dat = Array(comp_CH[para].data.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

		dat = Array(comp_CH[para].data.flat_top_cost[labelled_ind[train_ind], 2:end])
		mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

		mat = mat_c .* mat_ft

		# teind = [findall(i -> i == ("aug", test), labelled_ts[Not(train_ind)])[1] for test in tests]
	end

	gl1 = f[1:2, 1] = GridLayout()
	let
		axmain = Axis(gl1[1, 1], xlabel="training", ylabel="test")
		he = heatmap!(axmain, mat, colorscale=exp10, colormap=:oslo)

		cntmap = countmap(labelled_y[Not(train_ind)])

		lines!(axmain, [k[1]+0.5], [0.5, cntmap["C-LH"]+cntmap["C-EH"]+0.5], color=:white, linewidth=4)
		lines!(axmain, [k[1]+k[2]+0.5], [cntmap["C-LH"], length(labelled_y[Not(train_ind)])].+0.5, color=:white, linewidth=4)	

		lines!(axmain, [0, k[1]+k[2]].+0.5, [cntmap["C-LH"]+0.5], color=:white, linewidth=4)
		lines!(axmain, [k[1], length(labelled_y[train_ind])].+0.5, [cntmap["C-LH"]+cntmap["C-EH"]+0.5], color=:white, linewidth=4)
		
		Colorbar(gl1[1, 2], he, ticks=0:0.1:1)
	
		axmain.xticks = (1:length(labelled_ts[train_ind]), ["#$j" for (i, j) in labelled_ts[train_ind]])
		axmain.xticklabelrotation = π/2
		axmain.yticks = (1:length(labelled_ts[Not(train_ind)]), ["#$j" for (i, j) in labelled_ts[Not(train_ind)]])

		axmain.title = "C-LH 								C-EH 									CO"
	end
	
	gl2 = f[1:3, 2] = GridLayout()
	let
		axleftCL = Axis(gl2[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
			limits=(nothing, (0, 1.)), ylabel="C-LH", ylabelsize=15)
		axmidCL = Axis(gl2[1, 2], xticks=(0.7:1:0.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
			limits=(nothing, (0, 1.)))
		axrightCL = Axis(gl2[1, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)))

		ℓ = sum(train_ind)

		# CD = vcat([mat_c[i, i+1:end] for i in 1:ℓ-1]...)
		rng = 1:k[1]
		# println(rng)
		CD = vcat([mat_c[i, 1:end] for i in rng]...)
		hist!(axleftCL, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
		CD = vcat([mat_ft[i, 1:end] for i in rng]...)
		hist!(axmidCL, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
		CD = vcat([mat[i, 1:end] for i in rng]...)
		hist!(axrightCL, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
		if σ
			for (N, te) in enumerate(teind) 
				val, ind = findmax(mat[rng, te])
				println("C-LH train_ind: ", labelled_ts[train_ind][rng[1]+ind-1], ", K = $(round(val, digits=2))")
				hlines!(axleftCL, mat_c[rng[1]+ind-1, te], color=colors[N])
				hlines!(axmidCL, mat_ft[rng[1]+ind-1, te], color=colors[N])
				hlines!(axrightCL, val, color=colors[N])
			end
		end

		axleftCE = Axis(gl2[2, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
			limits=(nothing, (0, 1.)), ylabel="C-EH", ylabelsize=15)
		axmidCE = Axis(gl2[2, 2], xticks=(0.7:1:0.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
			limits=(nothing, (0, 1.)))
		axrightCE = Axis(gl2[2, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)))
		# CD = vcat([mat_ft[i, i+1:end] for i in 1:ℓ-1]...)
		# hist!(axmid, CD, normalization=:probability, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
		rng = (1+k[1]):(k[1]+k[2])
		# println(rng)
		CD = vcat([mat_c[i, 1:end] for i in rng]...)
		hist!(axleftCE, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray50)
		CD = vcat([mat_ft[i, 1:end] for i in rng]...)
		hist!(axmidCE, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray50)
		CD = vcat([mat[i, 1:end] for i in rng]...)
		hist!(axrightCE, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray50)
		if σ
			for (N, te) in enumerate(teind) 
				val, ind = findmax(mat[rng, te])
				println("C-EH train_ind: ", labelled_ts[train_ind][rng[1]+ind-1], ", K = $(round(val, digits=2))")
				println(mat_c[rng[1]+ind-1, te])
				hlines!(axleftCE, mat_c[rng[1]+ind-1, te], color=colors[N])
				hlines!(axmidCE, mat_ft[rng[1]+ind-1, te], color=colors[N])
				hlines!(axrightCE, val, color=colors[N])
			end
		end
		
		axleftCO = Axis(gl2[3, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
			limits=(nothing, (0, 1.)), ylabel="CO", ylabelsize=15)
		axmidCO = Axis(gl2[3, 2], xticks=(0.7:1:0.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
			limits=(nothing, (0, 1.)))
		axrightCO = Axis(gl2[3, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)))

		# CD = vcat([mat[i, i+1:end] for i in 1:ℓ-1]...)
		# hist!(axright, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
		rng = (1+k[1]+k[2]):(k[1]+k[2]+k[3])
		# println(rng)
		CD = vcat([mat_c[i, 1:end] for i in rng]...)
		hist!(axleftCO, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
		CD = vcat([mat_ft[i, 1:end] for i in rng]...)
		hist!(axmidCO, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
		CD = vcat([mat[i, 1:end] for i in rng]...)
		hist!(axrightCO, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)

		if σ
			for (N, te) in enumerate(teind) 
				val, ind = findmax(mat[rng, te])
				println("CO train_ind: ", labelled_ts[train_ind][rng[1]+ind-1], ", K = $(round(val, digits=2))")
				hlines!(axleftCO, mat_c[rng[1]+ind-1, te], color=colors[N])
				hlines!(axmidCO, mat_ft[rng[1]+ind-1, te], color=colors[N])
				hlines!(axrightCO, val, color=colors[N])
			end
		end
	end
	
	gl3 = f[3, 1] = GridLayout()
	let
		dat = Array(comp_CH[para].data.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

		dat = Array(comp_CH[para].data.flat_top_cost[labelled_ind[train_ind], 2:end])
		mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

		mat = mat_c .* mat_ft

		y = labelled_y[Not(train_ind)]
		ỹ, conf = svmpredict(model, mat)
		
		println("Accuracy: ", Accuracy()(ỹ, y), 
			"\nBalaced Accuracy: ", BalancedAccuracy()(ỹ, y), 
			"\nFalse Positive Rate: ", MulticlassFalsePositiveRate(;levels=["C-LH", "C-EH", "CO"], perm=[3,2,1])(ỹ, y))

		axmain = Axis3(gl3[1:2, 1:5], azimuth=-0.53π, xlabelvisible=false, ylabelvisible=false,
			zlabelvisible=false, viewmode=:stretch, elevation=0.05π)
		axleg = Axis(gl3[1:2, 0])
		hidedecorations!(axleg)
		hidespines!(axleg)

		scatter!(axmain, conf[1, :], conf[2, :], conf[3, :]; colormap=:tab10, 
			color=[col_dict[el] for el in ỹ],
			strokewidth = 0.1,
			label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
			markersize = [(i == j ? 20 : 10) for (i, j) in zip(y, ỹ)],
			marker = [(i == j ? '∘' : (:utriangle)) for (i, j) in zip(y, ỹ)]
		)	
		scatter!(axmain, [0], [0], [0], color=(:gray40, 0.3), markersize=20, strokewidth=0.2)
		
		plot_data = DataFrame(:ts => labelled_ts[Not(train_ind)],
			:predict => ỹ,
			:CH => y,
			:conf => [Tuple(el) for el in eachcol(conf)])
		
		pos = [(:right, :bottom), (:left, :bottom), (:right, :top), (:left, :top)]
		shuffle!(Random.seed!(3), plot_data)
		for (n, (ts, ỹ, y,	coord)) in enumerate(eachrow(plot_data))
			if y !== ỹ
				text!(axmain, coord..., text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=12, color=:gray40)
			elseif σ && in(ts[2], tests)
				text!(axmain, coord..., text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=12, color=:gray40)
			end	
		end

		group_marker = [MarkerElement(marker = marker, color = :black,
			strokecolor = :transparent, markersize = markersize) for (marker, markersize) in zip([:utriangle, '∘'], [12, 27])]

		Legend(gl3[1, 0], group_marker, ["incorrectly \n predicted", "correctly \n predicted"], 
			framevisible=false, labelsize=14, halign=:left, valign=:bottom)

		colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
		group_marker = [MarkerElement(marker = :circle, color = colors_3[i],
			strokecolor = :transparent, markersize = 12) for i in 1:3]

		Legend(gl3[2, 0], group_marker, ["CO", "C-EH", "C-LH"],
			framevisible=false, labelsize=14, halign=:left, valign=:top)

		axmain.xgridvisible = false
		axmain.ygridvisible = false
		axmain.zgridvisible = false
	end
	colsize!(gl3, 1, Relative(0.05))
	
	# if σ
		# gl4 = f[1, 3:4] = GridLayout()
		# let 
		# 	signal_1 = data_int.profile_data[("aug", train)]
		# 	signal_2 = data_int.profile_data[("aug", test)]
		# 	mat = dtw_cost_matrix(signal_2, signal_1, CosineDist(), transportcost=1.1)
		# 	cost, i1, i2 = DynamicAxisWarping.trackback(mat)
			
		# 	axmain1 = Axis(gl4[1:11, 3:14], title="cosine cost ($(round(cost, digits=1)))", titlesize=20)
		# 	begin
		# 		he = heatmap!(axmain1, mat, colormap=:thermal, colorscale=sqrt)
		# 		Colorbar(gl4[:, 15], he, ticks=0:5:100)

		# 		lines!(axmain1, i2, i1, color=:white, linewidth=3)
		# 		axmain1.xticklabelsvisible=false
		# 		axmain1.yticklabelsvisible=false
		# 		axmain1.xticks=1:13
		# 		axmain1.xticksvisible=false
		# 		axmain1.yticksvisible=false

		# 		for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
		# 			text!(axmain1, (i-0.1, j-0.7),
		# 				text="$(round(mat[i, j], digits=2))",
		# 				color=(:white, 0.5)
		# 			)
		# 		end
		# 		# for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
		# 		# 	text!(axmain1, (i-0.1, j-0.7),
		# 		# 		text="$(round(mat[i, j], digits=2))",
		# 		# 		color=(:white, 0.5)
		# 		# 	)
		# 		# end
		# 	end

		# 	axte = Axis(gl4[1:11, 1:2])
		# 	begin
		# 		for (n, feat) in enumerate(data_int.features)
		# 			lines!(axte, -signal_2[n, :][:]./feat_norm_dict[feat], 1:size(signal_2, 2), color=feat_color_dict[feat],
		# 				linestyle=feat_style_dict[feat])
		# 		end
		# 		axte.topspinevisible=false
		# 		axte.leftspinevisible=false
		# 		axte.rightspinevisible=false
		# 		axte.bottomspinevisible=false
		# 		axte.xgridvisible=false
		# 		axte.ygridvisible=false
		# 		axte.yticklabelsvisible=false
		# 		axte.xticklabelsvisible=false
		# 		axte.yticksvisible=false
		# 		axte.xticksvisible=false
		# 		axte.xlabel = "test (#$test)"
		# 		axte.xlabelsize = 17
		# 	end

		# 	axtr = Axis(gl4[12:14, 3:14])
		# 	begin
		# 		for (n, feat) in enumerate(data_int.features)
		# 			lines!(axtr, 1:size(signal_1, 2), signal_1[n, :]./feat_norm_dict[feat], color=feat_color_dict[feat],
		# 				linestyle=feat_style_dict[feat])
		# 		end
		# 		axtr.leftspinevisible=false
		# 		axtr.rightspinevisible=false
		# 		axtr.bottomspinevisible=false
		# 		axtr.xgridvisible=false
		# 		axtr.ygridvisible=false
		# 		axtr.yticklabelsvisible=false
		# 		axtr.xticklabelsvisible=false
		# 		axtr.yticksvisible=false
		# 		axtr.xticksvisible=false
		# 		axtr.ylabel = "train (#$train)"
		# 		axtr.ylabelsize = 17
		# 	end
		# end

		# gl5 = f[2, 3:4] = GridLayout()	
		# let 
		# 	signal_1 = data_int.flat_top_data[("aug", train)]
		# 	signal_2 = data_int.flat_top_data[("aug", test)]
		# 	mat = dtw_cost_matrix(signal_2, signal_1, Euclidean(), transportcost=1.1)
		# 	cost, i1, i2 = DynamicAxisWarping.trackback(mat)
			
		# 	axmain1 = Axis(gl5[1:11, 3:14], title="flat-top cost ($(round(cost, digits=1)))", titlesize=20)
		# 	begin
		# 		he = heatmap!(axmain1, mat, colormap=:thermal)
		# 		Colorbar(gl5[:, 15], he)

		# 		lines!(axmain1, i2, i1, color=:white, linewidth=3)
		# 		axmain1.xticklabelsvisible=false
		# 		axmain1.yticklabelsvisible=false
		# 		axmain1.xticks=1:13
		# 		axmain1.xticksvisible=false
		# 		axmain1.yticksvisible=false

		# 		for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
		# 			text!(axmain1, (i-0.1, j-0.7),
		# 				text="$(round(mat[i, j], digits=2))",
		# 				color=(:white, 0.5)
		# 			)
		# 		end
		# 		# for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
		# 		# 	text!(axmain1, (i-0.1, j-0.7),
		# 		# 		text="$(round(mat[i, j], digits=2))",
		# 		# 		color=(:white, 0.5)
		# 		# 	)
		# 		# end
		# 	end

		# 	axte = Axis(gl5[1:11, 1:2])
		# 	begin
		# 		for (n, feat) in enumerate(data_int.features)
		# 			lines!(axte, -signal_2[n, :][:]./feat_norm_dict[feat], 1:size(signal_2, 2), color=feat_color_dict[feat],
		# 				linestyle=feat_style_dict[feat])
		# 		end
		# 		axte.topspinevisible=false
		# 		axte.leftspinevisible=false
		# 		axte.rightspinevisible=false
		# 		axte.bottomspinevisible=false
		# 		axte.xgridvisible=false
		# 		axte.ygridvisible=false
		# 		axte.yticklabelsvisible=false
		# 		axte.xticklabelsvisible=false
		# 		axte.yticksvisible=false
		# 		axte.xticksvisible=false
		# 		axte.xlabel = "test (#$test)"
		# 		axte.xlabelsize = 17
		# 	end

		# 	axtr = Axis(gl5[12:14, 3:14])
		# 	begin
		# 		for (n, feat) in enumerate(data_int.features)
		# 			lines!(axtr, 1:size(signal_1, 2), signal_1[n, :]./feat_norm_dict[feat], color=feat_color_dict[feat],
		# 				linestyle=feat_style_dict[feat])
		# 		end
		# 		axtr.topspinevisible=false
		# 		axtr.leftspinevisible=false
		# 		axtr.rightspinevisible=false
		# 		axtr.bottomspinevisible=false
		# 		axtr.xgridvisible=false
		# 		axtr.ygridvisible=false
		# 		axtr.yticklabelsvisible=false
		# 		axtr.xticklabelsvisible=false
		# 		axtr.yticksvisible=false
		# 		axtr.xticksvisible=false
		# 		axtr.ylabel = "train (#$train)"
		# 		axtr.ylabelsize = 17
		# 	end
		# end
	# end
	# Label(f[0, :], text = "", fontsize = 20)
	# save("/Users/joe/Project/PhD/EuroFusion/EUROfusion_ML_2024/Presentations/Multidimensional_qunatity_classification_20_01_25/CO_EH_LH_classification.png", f)
	display(GLMakie.Screen(), f)
	# f
end

# plot the final results
let
	begin 
        trial = 4
        para = naming(["IP", "PNBI"])
        begin
            if trial == 1
                CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots

                COchoices = sort(tr_1[1:COn])
                CEchoices = sort(tr_1[COn+1:COn+CEn])
                CLchoices = sort(tr_1[COn+CEn+1:end])
            elseif trial == 2
                CLn, CEn, COn = trial_dataset_meta[:trial_2][2:4]
                tr_2 = (@subset trials_D @byrow :trial_2 == 1).shots

                COchoices = sort(tr_2[1:COn])
                CEchoices = sort(tr_2[COn+1:COn+CEn])
                CLchoices = sort(tr_2[COn+CEn+1:end])
            elseif trial == 3
                CLn, CEn, COn = trial_dataset_meta[:trial_3][2:4]
                tr_3 = (@subset trials_D @byrow :trial_3 == 1).shots

                COchoices = sort(tr_3[1:COn])
                CEchoices = sort(tr_3[COn+1:COn+CEn])
                CLchoices = sort(tr_3[COn+CEn+1:end])
            elseif trial == 4
                CLn, CEn, COn = trial_dataset_meta[:trial_4][2:4]
                tr_4 = (@subset trials_D @byrow :trial_4 == 1).shots

                COchoices = sort(tr_4[1:COn])
                CEchoices = sort(tr_4[COn+1:COn+CEn])
                CLchoices = sort(tr_4[COn+CEn+1:end])
            elseif trial == 5
                CLn, CEn, COn = trial_dataset_meta[:trial_5][2:4]
                tr_5 = (@subset trials_D @byrow :trial_5 == 1).shots

                COchoices = sort(tr_5[1:COn])
                CEchoices = sort(tr_5[COn+1:COn+CEn])
                CLchoices = sort(tr_5[COn+CEn+1:end])
            end
            labelled_data = let
                labelled_shots = vcat(CLchoices, CEchoices, COchoices)
                labels = vcat(repeat(["C-LH"], length(CLchoices)),
                    repeat(["C-EH"], length(CEchoices)),
                    repeat(["CO"], length(COchoices))
                ) 
                d = OrderedDict([ts for ts in labelled_shots] .=> labels)
            end
        end
        comp_CH = emp("complete_CH_trial_$(trial)_")
    end
	begin
		shots = [j for (i, j) in comp_CH[para].results.ts]
		labels = (@subset all_2D_features @byrow begin
			:tok == "aug"
			in(:shot, shots)
		end).current_heating
		comp_CH[para].results.conf = [col for col in eachcol(comp_CH[para].confidence)]
		curr_over = (@subset comp_CH[para].results @byrow :predict == "CO")
        
        EH_ind = findall(i -> in(i, ["C-EH"]), comp_CH[para].results.predict)
        LH_ind = findall(i -> in(i, ["C-LH"]), comp_CH[para].results.predict)
        CO_ind = findall(i -> in(i, ["CO"]), comp_CH[para].results.predict)

		acc_LH, acc_EH, acc_CO = Accuracy()(comp_CH[para].results.predict[LH_ind], labels[LH_ind]), Accuracy()(comp_CH[para].results.predict[EH_ind], labels[EH_ind]), Accuracy()(comp_CH[para].results.predict[CO_ind], labels[CO_ind]) 
        b_acc = BalancedAccuracy()(comp_CH[para].results.predict, labels)
		fpr = MulticlassFalsePositiveRate(;levels=["C-LH", "C-EH", "CO"], perm=[3,2,1])(comp_CH[para].results.predict, labels)
		println("accuracy LH = ($(acc_LH), $(acc_EH), $(acc_CO)), balanced accuracy = $(b_acc), false positive rate = $(fpr)")

		col_dict = Dict("CO"=>1, "C-EH"=>2, "C-LH"=>3)
		ỹ = comp_CH[para].results.predict
		y = labels

		f = Figure()
		a = Axis3(f[1,1])
	
		scatter!(comp_CH[para].confidence[1, :], comp_CH[para].confidence[2, :], comp_CH[para].confidence[3, :]; colormap=:tab10, 
			color=[col_dict[el] for el in ỹ],
			strokewidth = 0.1,
			label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
			markersize = [(i == j ? 20 : 10) for (i, j) in zip(y, ỹ)],
			marker = [(i == j ? '∘' : (:utriangle)) for (i, j) in zip(y, ỹ)],
            inspector_label=(self, i, pos)->"($(comp_CH[para].results.ts[i][2]), $(labels[i]))"
		)
        DataInspector(f)
		scatter!([0], [0], [0], color=(:gray40, 0.3), markersize=20, strokewidth=0.2)
		
		pos = [(:right, :bottom), (:left, :bottom), (:right, :top), (:left, :top)]
		Random.seed!(123)
		shuffle!(curr_over)
		for (n, (ts, coord)) in enumerate(zip(curr_over.ts, curr_over.conf))
			text!(a, coord..., text="$(ts[2])", align=pos[mod(n,4)+1], fontsize=10)
		end

		# extra = (@subset res @byrow in(:ts, [("aug", shot) for shot in [27930]]))
		# for (n, (ts, coord)) in enumerate(zip(extra.ts, extra.conf))
		#     text!(a, coord, text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=7)
		# end

		axislegend(a, unique=true, merge=true, position=:lb)
		a.xgridvisible = false
		a.ygridvisible = false
		a.zgridvisible = false
	end
	# display(GLMakie.Screen(), f)
	f
end
complete_CH_trial_1_["IP_PNBI"].results

# Determine train case used for classification
let 
	begin 
        trial = 4
        para = naming(["IP", "PNBI"])
        begin
            if trial == 1
                CLn, CEn, COn = trial_dataset_meta[:trial_1][2:4]
                tr_1 = (@subset trials_D @byrow :trial_1 == 1).shots

                COchoices = sort(tr_1[1:COn])
                CEchoices = sort(tr_1[COn+1:COn+CEn])
                CLchoices = sort(tr_1[COn+CEn+1:end])
            elseif trial == 2
                CLn, CEn, COn = trial_dataset_meta[:trial_2][2:4]
                tr_2 = (@subset trials_D @byrow :trial_2 == 1).shots

                COchoices = sort(tr_2[1:COn])
                CEchoices = sort(tr_2[COn+1:COn+CEn])
                CLchoices = sort(tr_2[COn+CEn+1:end])
            elseif trial == 3
                CLn, CEn, COn = trial_dataset_meta[:trial_3][2:4]
                tr_3 = (@subset trials_D @byrow :trial_3 == 1).shots

                COchoices = sort(tr_3[1:COn])
                CEchoices = sort(tr_3[COn+1:COn+CEn])
                CLchoices = sort(tr_3[COn+CEn+1:end])
            elseif trial == 4
                CLn, CEn, COn = trial_dataset_meta[:trial_4][2:4]
                tr_4 = (@subset trials_D @byrow :trial_4 == 1).shots

                COchoices = sort(tr_4[1:COn])
                CEchoices = sort(tr_4[COn+1:COn+CEn])
                CLchoices = sort(tr_4[COn+CEn+1:end])
            elseif trial == 5
                CLn, CEn, COn = trial_dataset_meta[:trial_5][2:4]
                tr_5 = (@subset trials_D @byrow :trial_5 == 1).shots

                COchoices = sort(tr_5[1:COn])
                CEchoices = sort(tr_5[COn+1:COn+CEn])
                CLchoices = sort(tr_5[COn+CEn+1:end])
            end
            labelled_data = let
                labelled_shots = vcat(CLchoices, CEchoices, COchoices)
                labels = vcat(repeat(["C-LH"], length(CLchoices)),
                    repeat(["C-EH"], length(CEchoices)),
                    repeat(["CO"], length(COchoices))
                ) 
                d = OrderedDict([ts for ts in labelled_shots] .=> labels)
            end
        end
        comp_CH = emp("complete_CH_trial_$(trial)_")
        shot = 39126
        best = 4
    end	
	begin
		HYP = comp_CH[para]

		data = DataFrame(HYP.kernel, ["$(ts)" for ts in HYP.results.ts])
		insertcols!(data, 1, :train => collect(keys(labelled_data)))

		ind = findfirst(i -> i == "(\"aug\", $(shot))", names(data))
		CL_train, CE_train, CO_train = data[1:CLn, ind], data[CLn+1:CLn+CEn, ind], data[CLn+CEn+1:end, ind]
		CL_dict = sort(Dict(CL_train .=> 1:CLn), byvalue=false, rev=true)
		CE_dict = sort(Dict(CE_train .=> 1:CEn), byvalue=false, rev=true)
		CO_dict = sort(Dict(CO_train .=> 1:COn), byvalue=false, rev=true)

		CL_ind = collect(values(CL_dict))[1:best]
		CE_ind = collect(values(CE_dict))[1:best] 
		CO_ind = collect(values(CO_dict))[1:best] 	

		data[vcat(CL_ind, CLn.+CE_ind, (CLn+CEn).+CO_ind), vcat(1, ind)]
	end
end

# See the closeness
let
    DTW_plot(("aug", 17389), ("aug", 26881), naming(["IP", "PNBI"]); section="IP")
end