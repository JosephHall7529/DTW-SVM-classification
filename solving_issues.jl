# common
labelled = (@subset trials_D @byrow :trial_1 == 1)
unseen = (@subset trials_D @byrow :trial_1 !== 1)
levels = ["LH", "EH", "CO"]
train, test = partition(labelled, 0.7; rng=12, stratify=labelled.label)

# old
data_dtw = let
    dict = OrderedDict()
    # collect(powerset(["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "LI", "NGW"], 1, 1))
    for para in [["IP"]]
        println(para)
        tss = tok_shots((@subset which(para) @byrow in(:tok, ["aug"])))
        dict[naming(para)] = DTW_hyp(tss, para, L=100) 
    end
    dict
end
data_dtw["IP"].flat_top_data[("aug", 42432)]

# hyperparameter optimise
## setup input data
begin
    cntmap=countmap(labelled.label)
    COchoices = (@subset labelled @byrow :label == "CO").shots
    CEchoices = (@subset labelled @byrow :label == "EH").shots
    CLchoices = (@subset labelled @byrow :label == "LH").shots
    data_dict_old = let
        labelled_shots = vcat(CLchoices, CEchoices, COchoices)
        labels = vcat(repeat(["LH"], length(CLchoices)),
            repeat(["EH"], length(CEchoices)),
            repeat(["CO"], length(COchoices))
        ) 
        d = OrderedDict([ts for ts in labelled_shots] .=> labels)
    end 
    k = [Int(div(cntmap["LH"], 10/7)), Int(div(cntmap["EH"], 10/7)), Int(div(cntmap["CO"], 10/7))]
    println(k)
end
# hyperparameter search
hyp_old = hyper_parameter_search(data_dtw["IP"], data_dict_old, k, interesting="CO", N=60, metric="MulticlassFalsePositiveRate(;levels=$levels, perm=[1, 2, 3], checks=false)", max=false)

# fit + predict
begin
    CC = 0.003
    FTC = 1.4
    
    tr_shots = train.shots
    train_names = df_ts_naming(tr_shots)
    tr_ind = [data_dtw["IP"].shot_dict[ts] for ts in tr_shots]
    X = exp.(-(Array(data_dtw["IP"].cosine_cost[tr_ind, train_names]) ./ CC).^2) .* 
        exp.(-(Array(data_dtw["IP"].flat_top_cost[tr_ind, train_names]) ./ FTC).^2)

    K_old = X 
    model = svmtrain(K_old, train.label, kernel=Kernel.Precomputed, cost=1.0)
    
    te_shots = test.shots
    test_names = df_ts_naming(te_shots)
    te_ind = [data_dtw["IP"].shot_dict[ts] for ts in te_shots]
    X = exp.(-(Array(data_dtw["IP"].cosine_cost[tr_ind, test_names]) ./ CC).^2) .* 
        exp.(-(Array(data_dtw["IP"].flat_top_cost[tr_ind, test_names]) ./ FTC).^2)
    KK_old = X

    ỹ_old, _ = svmpredict(model, KK_old)
    ỹ_old
    # acc_int[threadid()][i, j] += emp(metric)(ỹ, labelled_y[Not(train_ind)])
end

# testing
hyperparameters_meta!(hyp_old)
begin
    begin
        cntmap = countmap(values(hyp_old.labelled_data))
        CLn, CEn, COn = cntmap["LH"], cntmap["EH"], cntmap["CO"]
        labelled_ts = collect(hyp_old.labelled_data |> keys)
        labelled_y = collect(hyp_old.labelled_data |> values)
        shot_dict = Dict([a => n for (n, a) in enumerate(hyp_old.data.tok_shots)])
        labelled_ind = [shot_dict[k] for k in labelled_ts]
        labels = labelled_y |> unique

        reduce_ind = [shot_dict[ts] for ts in labelled_ts]
        reduced_name = df_ts_naming(labelled_ts)
        
        cos_, ft_, C_ = 0.003, 1.4, 1.0

        X = exp.(-(Array(hyp_old.data.cosine_cost[reduce_ind, reduced_name]) ./ (cos_)).^2) .* 
            exp.(-(Array(hyp_old.data.flat_top_cost[reduce_ind, reduced_name]) ./ (ft_)).^2) 
        
        ℓ_data = length(labelled_y)
        # k = [Int(div(*(CLn, 7), 10)), Int(div(*(CEn, 7), 10)), Int(div(*(COn, 7),10))]
    end

    Results = DataFrame()
    for S in 0:50
        begin
            train_ts, test_ts = partition(labelled_ts, 0.7, rng=S, stratify=coerce(labelled_y, OrderedFactor))
            train_ind = [in(i, train_ts) ? true : false for i in labelled_ts]

            model = svmtrain(X[train_ind, train_ind], labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C_)

            KK = X[train_ind, Not(train_ind)]
            ỹ, _ = svmpredict(model, KK)
            y = labelled_y[Not(train_ind)]
            # CO_EH_ind = findall(i -> in(i, ["CO", "C-EH"]), ỹ)
            
            acc = Accuracy()(ỹ, y)
            b_acc = BalancedAccuracy()(ỹ, y)

            fpr = MulticlassFalsePositiveRate(;levels=["LH", "EH", "CO"], perm=[1,2,3])(ỹ, y)
            # println("accuracy = $(acc), balanced accuracy = $(b_acc), false positive rate = $(fpr)")

            SVL, SVE, SVO = model.SVs.nSV
            # println(SVL, ", $(SVE), $(SVO)")
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
    Results
end






# new
# make optimisation structure fill database with profile/flattop data and create cosine/flattop cost
# model_new = DTW_SVM(["IP"]; target_shot_length=100, transportcost=1.1)
# fill_dtw_data!(model_new.database, 100)
# model_fill_train!(model_new, train.shots, train.label)
# model_fill_test!(model_new, test.shots)

# fit and predict
begin
    CC = 0.003
    FTC = 1.4

    train_fast = model_new.database.training_data
    tr_shots = train.shots
    train_names = df_ts_naming(tr_shots)
    ind = sortperm(train_fast.metadata, :order)
    X = exp.(-(Array(train_fast.cosine_cost[ind, train_names]) ./ CC).^2) .*
        exp.(-(Array(train_fast.flat_top_cost[ind, train_names]) ./ FTC).^2)

    K_new = X
    model = svmtrain(K_new, train.label, kernel=Kernel.Precomputed, cost=1.0)

    te_shots = test.shots
    test_names = df_ts_naming(te_shots)
    ind = sortperm(model_new.database.training_data.metadata, :order)

    X = exp.(-(Array(model_new.database.testing_data.cosine_cost[ind, test_names]) ./ CC).^2) .*
        exp.(-(Array(model_new.database.testing_data.flat_top_cost[ind, test_names]) ./ FTC).^2)
    
    KK_new = X

    ỹ_new, _ = svmpredict(model, KK_new)
    DataFrame(:pred_new => ỹ_new, :pred_old => ỹ_old, :actual => test.label)
end

# testing 
model_new = half_model(["IP"], labelled.shots, labelled.label;
    target_shot_length=500,
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
            tuning=MLJParticleSwarmOptimization.AdaptiveParticleSwarm(n_particles=10),
            # tuning=Grid(resolution=10),
            resampling=StratifiedCV(nfolds=4, shuffle=true),
            range=[r1, r2, r3],
            measure=MulticlassFalsePositiveRate(levels=levels),
            n=10,
            check_measure=true
    )

    mach = machine(self_tuning_tree, labelled.shots, labelled.label) 
    MLJBase.fit!(mach, verbosity=1)
end
mach_save = deepcopy(mach)
hyp_new = deepcopy(fitted_params(mach).best_model)

# hyp_new = deepcopy(fitted_params(mach).best_model)
# mach_new = deepcopy(mach)
begin
    Results = DataFrame()
    
    model_train = deepcopy(hyp_new.database.training_data)
    ℓ_data = size(model_train.metadata, 1)

    K = exp.(-(Array(model_train.cosine_cost[:, 2:end]) ./ (hyp_new.C_cosine)).^2) .* 
                exp.(-(Array(model_train.flat_top_cost[:, 2:end]) ./ (hyp_new.C_flat_top)).^2)
    for S in 0:100
        begin
            train_ind, test_ind = partition(1:ℓ_data, 0.3, rng=S, stratify=model_train.metadata.label, shuffle=true)
            sort!(train_ind)
            sort!(test_ind)

            model = svmtrain(K[train_ind, train_ind], model_train.metadata.label[train_ind], kernel=Kernel.Precomputed, cost=hyp_new.C_cost)

            ỹ, _ = svmpredict(model, K[train_ind, test_ind])
            y = model_train.metadata.label[test_ind]

            acc = Accuracy()(ỹ, y)
            b_acc = BalancedAccuracy()(ỹ, y)
            fpr = MulticlassFalsePositiveRate(;levels=["LH", "EH", "CO"])(ỹ, y)

            # println("accuracy = $(acc), balanced accuracy = $(b_acc), false positive rate = $(fpr)")

            SVL, SVE, SVO = model.SVs.nSV
            # println(SVL, ", $(SVE), $(SVO)")
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
    mean.(eachcol(Results[:, 2:4]))
    # hyp_old.hyperparameters_meta = Results
end
# hyperparameters_meta!(hyp_old)
eva = evaluate(hyp_new, [MulticlassFalsePositiveRate(levels=levels, perm=[1,2,3]), BalancedAccuracy(), Accuracy()])
hist(eva[:, 4], normalization=:pdf)