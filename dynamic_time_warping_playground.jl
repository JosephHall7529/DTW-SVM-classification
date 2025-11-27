using GLMakie
GLMakie.activate!()
colors = get(ColorSchemes.tab10, range(0,1,length=8));

trial = 5
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
# We now compare the results of each of the DTW structs specifically 
# CLn, CEn, COn = 20, 20, 10
# CLchoices = sort(sample(Random.seed!(13), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false))
# CEchoices = sort(sample(Random.seed!(13), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false))
# COchoices = sort(sample(Random.seed!(13), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false))

# complete_CH_trial_5 = Dict()
for trial in 1:5
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
	comp_CH = emp("complete_CH_trial_$trial")
	# if in(trial, [1,2,3])
	# 	params = [vcat(["IP"], [para]) for para in ["IP", "PECRH", "PICRH", "Q95", "LI", "NGW"]]
	# elseif trial == 4
	# 	params = [vcat(["NGW"], [para]) for para in ["IP", "PNBI", "PECRH", "PICRH", "Q95", "LI"]]
	# end
	for para in params
		para = naming(para)
		data_int = data_dtw[para]
		if !isempty(comp_CH)
			if in(para, collect(keys(comp_CH)))
				continue
			end
		end
		println("\n", para)
		k = [Int(div(CLn, 10/7)), Int(div(CEn, 10/7)), Int(div(COn, 10/7))]
		println(k)
		levels = ["C-LH", "C-EH", "CO"]
		
		comp_CH[para] = let
			labelled_data = let
				labelled_shots = vcat(CLchoices, CEchoices, COchoices)
				labels = vcat(repeat(["C-LH"], length(CLchoices)),
					repeat(["C-EH"], length(CEchoices)),
					repeat(["CO"], length(COchoices))
				) 
				d = OrderedDict([ts for ts in labelled_shots] .=> labels)
			end
			try 
				hyper_parameter_search(data_int, labelled_data, k, interesting="CO", N=60, metric="MulticlassFalsePositiveRate(;levels=$levels, perm=[1, 2, 3], checks=false)", max=false, cost=false)
			catch
				continue
			end
		end
	end
end

# Choosing C if desired
# for (para, hyp_struct) in complete_CH
# 	labelled_data = let
# 		labelled_shots = vcat(CLchoices, CEchoices, COchoices)
# 		labels = vcat(repeat(["C-LH"], length(CLchoices)),
# 			repeat(["C-EH"], length(CEchoices)),
# 			repeat(["CO"], length(COchoices))
# 		) 
# 		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
# 	end

# 	labelled_y = collect(labelled_data |> values)
# 	acc_int = [zeros(1) for _ in 1:nthreads()]
# 	labels = ["C-LH", "C-EH", "CO"]
# 	k = [6, 6, 6]

# 	Cs = [0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 1., 3., 5., 10., 30., 50., 100.]
        
# 	l = length(Cs)
# 	acc_int = [zeros(l, 1) for _ in 1:nthreads()]

# 	data = hyp_struct.data
# 	(cos_med, ft_med, C) = hyp_struct.hyperparameters

# 	Threads.@threads for S in ProgressBar(1:50)
# 		for i in 1:l
# 			C = Cs[i]
# 			train_ind = training_partion(labelled_data, labels, k, S=S)
# 			labelled_ind = [hyp_struct.data.shot_dict[k] for k in hyp_struct.labelled_ts]

# 			X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
# 				exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ ft_med).^2)

# 			K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
# 			model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C)

# 			KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
# 			ỹ, _ = svmpredict(model, KK)

# 			acc_int[threadid()][i] += MulticlassFalsePositiveRate(levels=["C-LH", "C-EH", "CO"], perm=[3,2,1])(ỹ, labelled_y[Not(train_ind)])
# 		end
# 	end
# 	Metric = map(+, acc_int...)
# 	Metric ./= 50
# 	hyp_struct.hyperparameters_meta.metric .= Metric
# end

# Table the results of the training data
begin
	trial = 1
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
	comp_CH = emp("complete_CH_trial_$trial")

	k = [Int(div(CLn, 10/7)), Int(div(CEn, 10/7)), Int(div(COn, 10/7))]
	train_tss = vcat(CLchoices, CEchoices, COchoices)
	n = length(train_tss)

	parameters = ["IP"]
	para = naming(parameters)
	data_int = data_dtw[para]

	labelled_data = let
		labelled_shots = vcat(CLchoices, CEchoices, COchoices)
		labels = vcat(repeat(["C-LH"], length(CLchoices)),
			repeat(["C-EH"], length(CEchoices)),
			repeat(["CO"], length(COchoices))) 
		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	labelled_ts = collect(labelled_data |> keys)
	labelled_y = collect(labelled_data |> values)
	# shot_dict = Dict([a => n for (n, a) in enumerate(data_int.tok_shots)])
	labelled_ind = [data_int.shot_dict[k] for k in labelled_ts]
	train_ind = training_partion(labelled_data, unique(labelled_y), k, S=S)

	X = exp.(-(Array(data_int.cosine_cost[!, 2:end]) ./ cos_med).^2) .*
		exp.(-(Array(data_int.flat_top_cost[!, 2:end]) ./ ft_med).^2)

	K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
	model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C)

	dat = Array(data_int.cosine_cost[labelled_ind[train_ind], 2:end])
	mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

	dat = Array(data_int.flat_top_cost[labelled_ind[train_ind], 2:end])
	mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

	mat = mat_c .* mat_ft

	teind = [findall(i -> i == ("aug", test), labelled_ts[Not(train_ind)])[1] for test in tests]
end

# Choose the model we want to plot the results of
para, S = "IP_PNBI", 123
println(para)
println(complete_CH[para].model.SVs.nSV)
# println(complete_CH[para].hyperparameters_meta.metric[1])
cos_med, ft_med, C = complete_CH[para].hyperparameters
data_int = complete_CH[para].data

labelled_data = let
	labelled_shots = vcat(CLchoices, CEchoices, COchoices)
	labels = vcat(repeat(["C-LH"], length(CLchoices)),
		repeat(["C-EH"], length(CEchoices)),
		repeat(["CO"], length(COchoices))) 
	d = OrderedDict([ts for ts in labelled_shots] .=> labels)
end
labelled_ts = collect(labelled_data |> keys)
labelled_y = collect(labelled_data |> values)

# Plot the train test case
let 
	# CairoMakie.activate!()
	f = Figure(size=(1000, 900))
	tests = [15714, 33371]
	opaque = ["IP", "PNBI"]
	σ = true 
	begin
		k = [Int(div(CLn, 10/7)), Int(div(CEn, 10/7)), Int(div(COn, 10/7))]
		train_tss = vcat(CLchoices, CEchoices, COchoices)
		n = length(train_tss)

		labelled_data = let
			labelled_shots = vcat(CLchoices, CEchoices, COchoices)
			labels = vcat(repeat(["C-LH"], length(CLchoices)),
				repeat(["C-EH"], length(CEchoices)),
				repeat(["CO"], length(COchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		# shot_dict = Dict([a => n for (n, a) in enumerate(data_int.tok_shots)])
		labelled_ind = [data_int.shot_dict[k] for k in labelled_ts]
		train_ind = training_partion(labelled_data, unique(labelled_y), k, S=S)

		X = exp.(-(Array(data_int.cosine_cost[!, 2:end]) ./ cos_med).^2) .*
			exp.(-(Array(data_int.flat_top_cost[!, 2:end]) ./ ft_med).^2)

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

		dat = Array(data_int.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

		dat = Array(data_int.flat_top_cost[labelled_ind[train_ind], 2:end])
		mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

		mat = mat_c .* mat_ft

		teind = [findall(i -> i == ("aug", test), labelled_ts[Not(train_ind)])[1] for test in tests]
	end

	gl1 = f[1:2, 1] = GridLayout()
	let
		axmain = Axis(gl1[1, 1])
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
		dat = Array(data_int.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

		dat = Array(data_int.flat_top_cost[labelled_ind[train_ind], 2:end])
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
	println(para)
	begin
		shots = [j for (i, j) in complete_CH[para].results.ts]
		labels = (@subset all_2D_features @byrow begin
			:tok == "aug"
			in(:shot, shots)
		end).current_heating
		complete_CH[para].results.conf = [col for col in eachcol(complete_CH[para].confidence)]
		curr_over = (@subset complete_CH[para].results @byrow :predict == "CO")

		acc, b_acc = Accuracy()(complete_CH[para].results.predict, labels), BalancedAccuracy()(complete_CH[para].results.predict, labels)
		fpr = MulticlassFalsePositiveRate(;levels=levels=["C-LH", "C-EH", "CO"], perm=[3,2,1])(complete_CH[para].results.predict, labels)
		println("accuracy = $(acc), balanced accuracy = $(b_acc), false positive rate = $(fpr)")

		col_dict = Dict("CO"=>1, "C-EH"=>2, "C-LH"=>3)
		ỹ = complete_CH[para].results.predict
		y = labels

		f = Figure()
		a = Axis3(f[1,1])
	
		scatter!(complete_CH[para].confidence[1, :], complete_CH[para].confidence[2, :], complete_CH[para].confidence[3, :]; colormap=:tab10, 
			color=[col_dict[el] for el in ỹ],
			strokewidth = 0.1,
			label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
			markersize = [(i == j ? 20 : 10) for (i, j) in zip(y, ỹ)],
			marker = [(i == j ? '∘' : (:utriangle)) for (i, j) in zip(y, ỹ)]
		)
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
# compare shots in all 8 parameters
let
	train = 37094
	test = 37093
	# 36177, 29772, 13622, 34610, 17282, 15714, 19314, 32135, 34441, 32464, 33407
	x1_lim = 10
	x2_lim = 10
	fig1 = Figure(size=(1000, 400))
	ax1 = Axis(fig1[1:3, 2:4], title="#$train", titlesize=15)
	ax2 = Axis(fig1[4:6, 2:4])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", train)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=param_colors, label=label)
	end
	ax1.limits = ((0., x1_lim), (-0.1, 2.5))
	ax2.limits = ((0., x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	linkxaxes!(ax1, ax2)
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)
	
	ax3 = Axis(fig1[1:3, 5:7], title="#$test", titlesize=15)
	ax4 = Axis(fig1[4:6, 5:7])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", test)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=param_colors, label=label)
	end
	ax3.limits = ((0., x2_lim), (-0.1, 2.5))
	ax4.limits = ((0., x2_lim), (-0.1, 2.5))
	ax3.yticks = 0:0.5:2.5
	ax4.yticks = 0:0.5:2.5
	linkxaxes!(ax3, ax4)
	# fig1[1:2, 8] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	# fig1[4:5, 8] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/shot_comparison.png", fig1)
	display(GLMakie.Screen(), fig1)
end

results = let 
	k = [6, 6, 6]
	train_tss = vcat(CLchoices, CEchoices, COchoices)
	n = length(train_tss)

	labelled_data = let
		labelled_shots = vcat(CLchoices, CEchoices, COchoices)
		labels = vcat(repeat(["C-LH"], length(CLchoices)),
			repeat(["C-EH"], length(CEchoices)),
			repeat(["CO"], length(COchoices))) 
		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	labelled_ts = collect(labelled_data |> keys)
	labelled_y = collect(labelled_data |> values)

	next = DataFrame()
	for (para, hyp_struct) in complete_CH
		cos_med, ft_med, C = hyp_struct.hyperparameters
		parameters = Vector{String}()
		Cs = Vector{Number}()
		cos_hists = Vector{Tuple{Float64, Float64}}()
		cos_modes = Vector{Float64}()
		ft_hists = Vector{Tuple{Float64, Float64}}()
		ft_modes = Vector{Float64}()
		kernel_hists = Vector{Tuple{Float64, Float64}}()
		kernel_modes = Vector{Float64}()

		data_int = hyp_struct.data
		shot_dict = hyp_struct.data.shot_dict
		res = hyp_struct.results
		confidence = hyp_struct.confidence
		kernel_matrix = hyp_struct.kernel
		model = hyp_struct.model
		
		labelled_ind = [shot_dict[k] for k in labelled_ts]
		for S in 1:50
			train_ind = training_partion(labelled_data, unique(labelled_y), k, S=S)
			ℓ = sum(train_ind)

			X = exp.(-(Array(data_int.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
				exp.(-(Array(data_int.flat_top_cost[!, 2:end]) ./ ft_med).^2
			)

			K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
			model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C)

			dat = Array(data_int.cosine_cost[labelled_ind[train_ind], 2:end])
			mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

			CD = vcat([mat_c[i, i+1:end] for i in 1:ℓ-1]...)
			push!(cos_hists, extrema(CD))
			push!(cos_modes, mode(CD))

			dat = Array(data_int.flat_top_cost[labelled_ind[train_ind], 2:end])
			mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

			CD = vcat([mat_ft[i, i+1:end] for i in 1:ℓ-1]...)
			push!(ft_hists, extrema(CD))
			push!(ft_modes, mode(CD))

			mat = mat_c .* mat_ft

			CD = vcat([mat[i, i+1:end] for i in 1:ℓ-1]...)
			push!(kernel_hists, extrema(CD))
			push!(kernel_modes, mode(CD))
			# push!(parameters, hyp_search.data.parameter)
			# push!(Cs, C)
		end
		append!(next, DataFrame(:parameters => para,
			:C => C,
			:cos_hist_lb => [i for (i, j) in cos_hists],
			:cos_hist_ub => [j for (i, j) in cos_hists],
			:cos_mode => cos_modes,
			:ft_hist_lb => [i for (i, j) in ft_hists],
			:ft_hist_ub => [j for (i, j) in ft_hists],
			:ft_mode => ft_modes,
			:kernel_hist_lb => [i for (i, j) in kernel_hists],
			:kernel_hist_ub => [j for (i, j) in kernel_hists],
			:kernel_mode => kernel_modes))
	end
	next
end
results
# next = let 
# 	if isempty(results)
# 		return "results empty"
# 	end
# 	Bchoices = [("aug", shot) for shot in chosen_baseline_training]	
# 	k = [10, 10]
# 	train_tss = vcat(Bchoices, Hchoices)
# 	n = length(train_tss)

# 	labelled_data = let
# 		labelled_shots = vcat(Bchoices, Hchoices)
# 		labels = vcat(repeat(["baseline"], length(Bchoices)),
# 			repeat(["hybrid"], length(Hchoices))) 
# 		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
# 	end
# 	labelled_ts = collect(labelled_data |> keys)
# 	labelled_y = collect(labelled_data |> values)
# 	shot_dict = Dict([a => n for (n, a) in enumerate(data_int.tok_shots)])
# 	labelled_ind = [shot_dict[k] for k in labelled_ts]

# 	next = DataFrame()
# 	for (para, cos_med, ft_med) in eachrow(chosen_meta_FP[!, [:parameters, :COS, :FT]])
# 		if in(para, unique(results.parameters))
# 			continue
# 		end
# 		cos_hists = Vector{Tuple{Float64, Float64}}()
# 		cos_modes = Vector{Float64}()	
# 		ft_hists = Vector{Tuple{Float64, Float64}}()
# 		ft_modes = Vector{Float64}()	
# 		kernel_hists = Vector{Tuple{Float64, Float64}}()
# 		kernel_modes = Vector{Float64}()	
# 		for S in 1:100
# 			train_ind = training_partion(labelled_data, unique(labelled_y), k, S=S)
# 			ℓ = sum(train_ind)

# 			data_int = data_bhind[para]
# 			# data_int = data_bhind[*(sort(["PICRH", "PNBI"]), delim="_")]
# 			res, confidence, kernel_matrix, model = let
# 				labelled_data = let
# 					labelled_shots = vcat(Bchoices, Hchoices)
# 					labels = vcat(repeat(["baseline"], length(Bchoices)),
# 						repeat(["hybrid"], length(Hchoices))) 
# 					d = OrderedDict([ts for ts in labelled_shots] .=> labels)
# 				end
# 				classify!(data_int, labelled_data, (cos_med, ft_med); interesting="hybrid")
# 			end

# 			X = exp.(-(Array(data_int.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
# 				exp.(-(Array(data_int.flat_top_cost[!, 2:end]) ./ ft_med).^2)

# 			K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
# 			model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed)

# 			dat = Array(data_int.cosine_cost[labelled_ind[train_ind], 2:end])
# 			mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

# 			CD = vcat([mat_c[i, i+1:end] for i in 1:ℓ-1]...)
# 			push!(cos_hists, extrema(CD))
# 			push!(cos_modes, mode(CD))

# 			dat = Array(data_int.flat_top_cost[labelled_ind[train_ind], 2:end])
# 			mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

# 			CD = vcat([mat_ft[i, i+1:end] for i in 1:ℓ-1]...)
# 			push!(ft_hists, extrema(CD))
# 			push!(ft_modes, mode(CD))

# 			mat = mat_c .* mat_ft

# 			CD = vcat([mat[i, i+1:end] for i in 1:ℓ-1]...)
# 			push!(kernel_hists, extrema(CD))
# 			push!(kernel_modes, mode(CD))
# 		end
# 		append!(next, DataFrame(:parameters => para,
# 			:cos_hist_lb => [i for (i, j) in cos_hists],
# 			:cos_hist_ub => [j for (i, j) in cos_hists],
# 			:cos_mode => cos_modes,
# 			:ft_hist_lb => [i for (i, j) in ft_hists],
# 			:ft_hist_ub => [j for (i, j) in ft_hists],
# 			:ft_mode => ft_modes,
# 			:kernel_hist_lb => [i for (i, j) in kernel_hists],
# 			:kernel_hist_ub => [j for (i, j) in kernel_hists],
# 			:kernel_mode => kernel_modes))
# 	end
# 	next
# end
# append!(results, next) 

complete_CH["PNBI"].data.tok_shots
# SD = complete_CH["PNBI"].data.shot_dict
# complete_CH["PNBI"].labelled_ts
# yts = complete_CH["PNBI"].data.tok_shots[Not([SD[k] for k in complete_CH["PNBI"].labelled_ts])]
# [current_heating(ts...) for ts in yts]
complete_CH["PNBI"].hyperparameters_meta
N = 10
begin
	K = 50
	para = unique(results.parameters)[N]
	data = complete_CH[para]

	plot_D = (@subset results @byrow :parameters == para)
	fig = Figure();
	axleft = Axis(fig[1, 1], xticks=(1:0.2:1, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
		limits=(nothing, (0, 1.)), title="σ=$(round(complete_CH[para].hyperparameters[1], digits=1))")
	axmid = Axis(fig[1, 2], xticks=(1:0.2:1, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
		limits=(nothing, (0, 1.)), title="σ=$(round(complete_CH[para].hyperparameters[2], digits=1))")
	axright = Axis(fig[1, 3], xticks=(1:1:1, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)),
		title="σ=$(round(complete_CH[para].hyperparameters[3], digits=1))")
	
	boxplot!(axleft, ones(K), plot_D.cos_hist_lb, color=:gray70)
	boxplot!(axleft, ones(K), plot_D.cos_hist_ub, color=:gray50)
	boxplot!(axleft, ones(K), plot_D.cos_mode, color=:red)

	boxplot!(axmid, ones(K), plot_D.ft_hist_lb, color=:gray70)
	boxplot!(axmid, ones(K), plot_D.ft_hist_ub, color=:gray50)
	boxplot!(axmid, ones(K), plot_D.ft_mode, color=:red)

	boxplot!(axright, ones(K), plot_D.kernel_hist_lb, color=:gray70)
	boxplot!(axright, ones(K), plot_D.kernel_hist_ub, color=:gray50)
	boxplot!(axright, ones(K), plot_D.kernel_mode, color=:red)

	println(N, " ", para, " ", Int.(data.model.SVs.nSV))
	SD = data.data.shot_dict
	yts = data.data.tok_shots[Not([SD[k] for k in data.labelled_ts])]
	D = DataFrame(:pred => complete_CH[para].results.predict, :label => [current_heating(ts...) for ts in yts])
	fpr = (@subset data.hyperparameters_meta @byrow :C == C).metric[1]

	println("BA: ", round(BalancedAccuracy()(D.pred, D.label), digits=2), " FPR: ", round(fpr, digits=2))
	CO, CE, CL = (@subset D @byrow :label == "CO"), (@subset D @byrow :label == "C-EH"), (@subset D @byrow :label == "C-LH")
	println("CO accuracy: ", Accuracy()(CO.pred, CO.label))
	println("C-EH accuracy: ", Accuracy()(CE.pred, CE.label))
	println("C-LH accuracy: ", Accuracy()(CL.pred, CL.label))
	N += 1
	fig
end





















# HYBRID CLASSIFICATION

# We now compare the results of each of the DTW structs specifically 

chosen_baseline_training = [6483, 7558, 7643, 8007, 8093, 10729, 11126, 13603, 15714, 17389, 20280,  
	25719, 25757, 29623, 32225, 33303, 34489, 34490, 35975, 40203, 40411
]
Hchoices = sort([("aug", shot) for shot in [11190, 16688, 16736, 18046, 18869, 18880, 19314, 25764, 26338, 26913, 27930, 34769, 34770, 
	34774, 34775, 34776, 36443, 36408]])

Hn = length(Hchoices)

Bchoices = [("aug", shot) for shot in chosen_baseline_training]
Bn = length(Bchoices)

# complete_BH = Dict()
for para in BH_paras
	data = data_dtw[naming(para)]
	k = [10, 10]
	
	if !isempty(complete_BH)
		if in(naming(para), collect(keys(complete_BH)))
			continue
		end
	end

	complete_BH[naming(para)] = let
		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		try
			levels = ["baseline", "hybrid"]
			hyper_parameter_search(data, labelled_data, k, interesting="hybrid", N=50, metric="FalsePositiveRate(;levels=$levels)", max=false)
		catch
			continue
		end
	end

	if !in(naming(para), collect(keys(complete_BH)))
		continue
	end
	RES = complete_BH[naming(para)]
	cos_med, ft_med, C = complete_BH[naming(para)].hyperparameters

	train_tss = vcat(Bchoices, Hchoices)
	n = length(train_tss)
	ind = [data.shot_dict[ts] for ts in train_tss]

	dat = data.cosine_cost
	mat = exp.(-(Array(dat[ind, ind.+1]) ./ cos_med).^2)
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

	cos_mode = mode(CD)
	cos_ext = extrema(CD)

	dat = data.flat_top_cost
	mat = exp.(-(Array(dat[ind, ind.+1]) ./ ft_med).^2)
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

	ft_mode = mode(CD)
	ft_ext = extrema(CD)

	shots = [j for (i, j) in RES.results.ts]
	labels = (@subset all_2D_features @byrow begin
		:tok == "aug"
		in(:shot, shots)
	end).current_heating

	BACC = 0
	ACC_hyb = 0
	ACC_base = 0
	for S in 1:50	
		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
		labelled_ind = [shot_dict[k] for k in labelled_ts]

		shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
		labels = collect(labelled_data |> values) |> unique

		train_ind = training_partion(labelled_data, labels, k, S=S)

		X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
			exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ ft_med).^2)

		K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
		model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed, cost=C)

		KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
		ỹ, _ = svmpredict(model, KK)
		
		int_D = DataFrame(:pred => ỹ, :label => labelled_y[Not(train_ind)])

		BACC += BalancedAccuracy()(int_D.pred, int_D.label)
		hyb_D = (@subset int_D @byrow :label == "hybrid")
		ACC_hyb += Accuracy()(hyb_D.pred, hyb_D.label)
		base_D = (@subset int_D @byrow :label == "baseline")
		ACC_base += Accuracy()(base_D.pred, base_D.label)
	end
	BACC = BACC / 50
	ACC_hyb = ACC_hyb / 50
	ACC_base = ACC_base /50

	append!(RES.hyperparameters_meta, 
		DataFrame(:parameters => *(sort(para), delim="_"),
			:metric => RES.acc,
			:Bacc => BACC,
			:ACC_hyb => ACC_hyb,
			:ACC_base => ACC_base,
			:COS => cos_med,
			:cos_hist_mode => cos_mode,
			:cos_hist_ext => "$cos_ext",
			:FT => ft_med,
			:ft_hist_mode => ft_mode,
			:ft_hist_ext => "$ft_ext",
			:C => C,
			:nSV => "$(Int.(RES.model.SVs.nSV))"
			))
end
D = DataFrame()
for para in BH_paras

	if !in(naming(para), collect(keys(complete_BH)))
		continue
	end
	append!(D, complete_BH[naming(para)].hyperparameters_meta)
end
D |> vscodedisplay
# chosen_meta_FP = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/meta/chosen_baseline_false_positive_results.csv", DataFrame, stringtype=String)
# chosen_meta_BA = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/meta/chosen_baseline_balanced_accuracy_results.csv", DataFrame, stringtype=String)

# Choose the model we want to plot the results of
S = 1
n = 14
para = D[n, :parameters]
Bchoices = [("aug", shot) for shot in chosen_baseline_training]	
cos_med, ft_med, C = complete_BH[para].hyperparameters

data_int = data_dtw[para]
res, confidence, kernel_matrix, model = let
	labelled_data = let
		labelled_shots = vcat(Bchoices, Hchoices)
		labels = vcat(repeat(["baseline"], length(Bchoices)),
			repeat(["hybrid"], length(Hchoices))) 
		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	classify!(data_int, labelled_data, (cos_med, ft_med, C); interesting="hybrid")
end

# Plot the train test case
train = 8041
# test = 34490
test = 8042
let 
	GLMakie.activate!()
	F = Figure(size=(1600, 1150), title=para)
	opaque = ["NGW", "BETAPOL"]
	σ = false

	f = F[1, 1] = GridLayout()

	begin
		k = [10, 10]
		train_tss = vcat(Bchoices, Hchoices)
		n = length(train_tss)

		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		shot_dict = Dict([a => n for (n, a) in enumerate(data_int.tok_shots)])
		labelled_ind = [shot_dict[k] for k in labelled_ts]
		train_ind = training_partion(labelled_data, unique(labelled_y), k, S=S)

		X = exp.(-(Array(data_int.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
			exp.(-(Array(data_int.flat_top_cost[!, 2:end]) ./ ft_med).^2)

		K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
		model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed)

		col_dict = Dict("baseline"=>1, "hybrid"=>2, "training"=>3)
		CO_dict = Dict("CO"=> :cross, "C-EH"=> :ltriangle, "C-LH"=> :rtriangle)
		feat_norm_dict = Dict("IP" => 1, "PNBI" => 10, "PECRH" => 10, "PICRH" => 10, 
			"Q95" => 10, "NGW" => 1, "LI" => 1, "BETAPOL" => 1)
		feat_color_dict = Dict(
			"IP" => (colors[1], in("IP", opaque) ? 1 : 0.3), 
			"PNBI" => (colors[2], in("PNBI", opaque) ? 1 : 0.3), 
			"PECRH" => (colors[3], in("PECRH", opaque) ? 1 : 0.3), 
			"PICRH" => (colors[4], in("PICRH", opaque) ? 1 : 0.3), 
			"Q95" => (colors[5], in("Q95", opaque) ? 1 : 0.3), 
			"NGW" => (colors[6], in("NGW", opaque) ? 1 : 0.3), 
			"LI" => (colors[7], in("LI", opaque) ? 1 : 0.3), 
			"BETAPOL" => (colors[8], in("BETAPOL", opaque) ? 1 : 0.3))
		feat_style_dict = Dict(
			"IP" => :solid, 
			"PNBI" => :dash, 
			"PECRH" => :dash, 
			"PICRH" => :dash, 
			"Q95" => :dot, 
			"NGW" => :dot, 
			"LI" => :dash, 
			"BETAPOL" => :dash)

		dat = Array(data_int.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

		dat = Array(data_int.flat_top_cost[labelled_ind[train_ind], 2:end])
		mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

		mat = mat_c .* mat_ft

		if σ
			trind = findall(i -> i == ("aug", train), labelled_ts[train_ind])[1]
			teind = findall(i -> i == ("aug", test), labelled_ts[Not(train_ind)])[1]
		end
	end

	gl1 = f[1, 1] = GridLayout()
	let
		axmain = Axis(gl1[1, 1])
		he = heatmap!(axmain, mat, colorscale=exp10, colormap=:oslo)

		cntmap = countmap(labelled_y[Not(train_ind)])

		lines!(axmain, [k[1]+0.5], [0.5, cntmap["baseline"]+cntmap["hybrid"]+0.5], color=:white, linewidth=4)

		lines!(axmain, [0, k[1]+k[2]].+0.5, [cntmap["baseline"]+0.5], color=:white, linewidth=4)
		
		Colorbar(gl1[1, 2], he, ticks=0:0.1:1)
	
		axmain.xticks = (1:length(labelled_ts[train_ind]), ["#$j" for (i, j) in labelled_ts[train_ind]])
		axmain.xticklabelsize = 17
		axmain.yticklabelsize = 17
		axmain.xticklabelrotation = π/2
		axmain.yticks = (1:length(labelled_ts[Not(train_ind)]), ["#$j" for (i, j) in labelled_ts[Not(train_ind)]])

		axmain.title = "baseline                                                          hybrid"
		axmain.titlesize = 20
	end

	gl2 = f[1, 2] = GridLayout()
	let
		axleft = Axis(gl2[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
			limits=(nothing, (0, 1.)), title="σ=$(round(cos_med, digits=1))")
		axmid = Axis(gl2[1, 2], xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
			limits=(nothing, (0, 1.)), title="σ=$(round(ft_med, digits=1))")
		axright = Axis(gl2[1, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)))

		ℓ = sum(train_ind)

		CD = vcat([mat_c[i, i+1:end] for i in 1:ℓ-1]...)
		hist!(axleft, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
		if σ
			hlines!(axleft, [mat_c[trind, teind]], color=:red)
		end

		CD = vcat([mat_ft[i, i+1:end] for i in 1:ℓ-1]...)
		hist!(axmid, CD, normalization=:probability, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
		if σ
			hlines!(axmid, [mat_ft[trind, teind]], color=:red)
		end
		
		mat = mat_c .* mat_ft
		CD = vcat([mat[i, i+1:end] for i in 1:ℓ-1]...)
		hist!(axright, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
		if σ
			hlines!(axright, [mat[trind, teind]], color=:red)
		end
	end

	gl3 = f[2, 1:2] = GridLayout()
	let
		dat = Array(data_int.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

		dat = Array(data_int.flat_top_cost[labelled_ind[train_ind], 2:end])
		mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

		mat = mat_c .* mat_ft

		y = labelled_y[Not(train_ind)]
		ỹ, conf = svmpredict(model, mat)

		println("Accuracy: ", Accuracy()(ỹ, y), "\nBalaced Accuracy: ", BalancedAccuracy()(ỹ, y))

		ax_base_tr = Axis(gl3[1, 2:3], xticklabelsvisible=false, yticklabelsvisible=false, limits=(nothing, (-1, 1)))
		hidespines!(ax_base_tr)
		hidedecorations!(ax_base_tr)
		axmain = Axis(gl3[2:3, 2:3], xlabelvisible=false, ylabelvisible=false, limits=(nothing, (-2, 2)))
		axmain.xticks=(10000:10000:40000, ["10,000", "20,000", "30,000", "40,000"])
		ax_hyb_tr = Axis(gl3[4, 2:3], xticklabelsvisible=false, yticklabelsvisible=false, limits=(nothing, (-1, 1)))
		hidespines!(ax_hyb_tr)
		hidedecorations!(ax_hyb_tr)

		axleg = Axis(gl3[1:2, 1])
		hidedecorations!(axleg)
		hidespines!(axleg)

		te_x = [j for (i,j) in labelled_ts[Not(train_ind)]]
		hlines!(axmain, [0], color=:gray80, linestyle=:dash)
		scatter!(axmain, te_x, conf[1, :]; colormap=:tab10, 
			color=[col_dict[el] for el in ỹ],
			colorrange=(1, 3),
			strokewidth = 0.1,
			label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
			markersize = [(i == j ? 25 : 15) for (i, j) in zip(y, ỹ)],
			marker = [(i == j ? '∘' : (:utriangle)) for (i, j) in zip(y, ỹ)]
		)	
		
		plot_data = DataFrame(:ts => labelled_ts[Not(train_ind)],
			:predict => ỹ,
			:CH => y,
			:conf => [Tuple(el) for el in eachcol(conf)]
		)

		training_D = DataFrame(:ts => labelled_ts[train_ind],
			:label => labelled_y[train_ind],
			:conf => vcat([(ts[2], (0.6*cos(n+1))) for (n, ts) in enumerate(labelled_ts[train_ind])])
		)
		training_D.CH = [current_heating(ts...) for ts in training_D.ts]
		tr_base, tr_hyb = (@subset training_D @byrow :label == "baseline"), (@subset training_D @byrow :label == "hybrid")

		scatter!(ax_base_tr, tr_base.conf; colormap=:tab10, colorrange=(1,3), color=3, marker=[CO_dict[el] for el in tr_base.CH], markersize=17)
		scatter!(ax_hyb_tr, tr_hyb.conf; colormap=:tab10, colorrange=(1,3), color=3, marker=[CO_dict[el] for el in tr_hyb.CH], markersize=17)

		pos = shuffle(Random.seed!(2), [(:left, :bottom), (:right, :bottom), (:right, :top), (:left, :top)])
		shuffle!(Random.seed!(13), plot_data)
		for (n, (ts, coord)) in enumerate(zip(plot_data.ts, plot_data.conf))
			text!(axmain, ts[2], coord[1], text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=18, color=:gray40)
		end
		for (n, (ts, coord)) in enumerate(zip(tr_base.ts, tr_base.conf))
			text!(ax_base_tr, coord..., text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=18, color=:gray40)
		end
		for (n, (ts, coord)) in enumerate(zip(tr_hyb.ts, tr_hyb.conf))
			text!(ax_hyb_tr, coord..., text="$(ts[2])", align=pos[mod(n+4,4)+1], fontsize=18, color=:gray40)
		end

		begin
			group_marker = [MarkerElement(marker = marker, color = :black, 
				strokecolor = :transparent, markersize = markersize) for (marker, markersize) in zip([:utriangle, '∘'], [12, 27])
			]
			CH_marker = [MarkerElement(marker = marker, color = :black,
				strokecolor = :transparent, markersize = 12) for marker in [:ltriangle, :rtriangle, :cross]]
			colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3)
			)
			group_color = [MarkerElement(marker = :circle, color = colors_3[i],
				strokecolor = :transparent, markersize = 12) for i in 1:3
			]
			Legend(gl3[2:3, 1], [group_marker, group_color, CH_marker],
				[["incorrectly \npredicted", "correctly \npredicted"], ["baseline", "hybrid", "training"], ["C-EH", "C-LH", "CO"]], 
				["", "", ""], framevisible=false, labelsize=20, tellwidth=true, tellheight=true
			)
		end

		axmain.xgridvisible = false
		axmain.ygridvisible = false
		linkxaxes!(ax_base_tr, ax_hyb_tr, axmain)
	end

	if σ
		gl4 = f[1, 3:4] = GridLayout()
		let 
			signal_1 = data_int.profile_data[("aug", train)]
			signal_2 = data_int.profile_data[("aug", test)]
			mat = dtw_cost_matrix(signal_2, signal_1, CosineDist(), transportcost=1.1)
			cost, i1, i2 = DynamicAxisWarping.trackback(mat)
			
			axmain1 = Axis(gl4[1:11, 3:14], title="cosine cost ($(round(cost, digits=1)))  | σ = $(round(cos_med, digits=1))", titlesize=20)
			begin
				he = heatmap!(axmain1, mat, colormap=:thermal, colorscale=sqrt)
				Colorbar(gl4[:, 15], he, ticks=0:5:100)

				lines!(axmain1, i2, i1, color=:white, linewidth=3)
				axmain1.xticklabelsvisible=false
				axmain1.yticklabelsvisible=false
				axmain1.xticks=1:13
				axmain1.xticksvisible=false
				axmain1.yticksvisible=false

				for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
					text!(axmain1, (i-0.1, j-0.7),
						text="$(round(mat[i, j], digits=2))",
						color=(:white, 0.5)
					)
				end
				# for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
				# 	text!(axmain1, (i-0.1, j-0.7),
				# 		text="$(round(mat[i, j], digits=2))",
				# 		color=(:white, 0.5)
				# 	)
				# end
			end

			axte = Axis(gl4[1:11, 1:2])
			begin
				for (n, feat) in enumerate(data_int.features)
					lines!(axte, -signal_2[n, :][:]./feat_norm_dict[feat], 1:size(signal_2, 2), color=feat_color_dict[feat],
						linestyle=feat_style_dict[feat])
				end
				axte.topspinevisible=false
				axte.leftspinevisible=false
				axte.rightspinevisible=false
				axte.bottomspinevisible=false
				axte.xgridvisible=false
				axte.ygridvisible=false
				axte.yticklabelsvisible=false
				axte.xticklabelsvisible=false
				axte.yticksvisible=false
				axte.xticksvisible=false
				axte.xlabel = "test (#$test)"
				axte.xlabelsize = 17
			end

			axtr = Axis(gl4[12:14, 3:14])
			begin
				for (n, feat) in enumerate(data_int.features)
					lines!(axtr, 1:size(signal_1, 2), signal_1[n, :]./feat_norm_dict[feat], color=feat_color_dict[feat],
						linestyle=feat_style_dict[feat])
				end
				axtr.topspinevisible=false
				axtr.leftspinevisible=false
				axtr.rightspinevisible=false
				axtr.bottomspinevisible=false
				axtr.xgridvisible=false
				axtr.ygridvisible=false
				axtr.yticklabelsvisible=false
				axtr.xticklabelsvisible=false
				axtr.yticksvisible=false
				axtr.xticksvisible=false
				axtr.ylabel = "train (#$train)"
				axtr.ylabelsize = 17
			end
		end

		gl5 = f[2, 3:4] = GridLayout()	
		let 
			signal_1 = data_int.flat_top_data[("aug", train)]
			signal_2 = data_int.flat_top_data[("aug", test)]
			mat = dtw_cost_matrix(signal_2, signal_1, Euclidean(), transportcost=1.1)
			cost, i1, i2 = DynamicAxisWarping.trackback(mat)
			
			axmain1 = Axis(gl5[1:11, 3:14], title="flat-top cost ($(round(cost, digits=1)))  | σ = $(round(ft_med, digits=1))", titlesize=20)
			begin
				he = heatmap!(axmain1, mat, colormap=:thermal)
				Colorbar(gl5[:, 15], he)

				lines!(axmain1, i2, i1, color=:white, linewidth=3)
				axmain1.xticklabelsvisible=false
				axmain1.yticklabelsvisible=false
				axmain1.xticks=1:13
				axmain1.xticksvisible=false
				axmain1.yticksvisible=false

				for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
					text!(axmain1, (i-0.1, j-0.7),
						text="$(round(mat[i, j], digits=2))",
						color=(:white, 0.5)
					)
				end
				# for (i, j) in Iterators.product(2:13:size(signal_1, 2), 1:13:size(signal_2, 2))
				# 	text!(axmain1, (i-0.1, j-0.7),
				# 		text="$(round(mat[i, j], digits=2))",
				# 		color=(:white, 0.5)
				# 	)
				# end
			end

			axte = Axis(gl5[1:11, 1:2])
			begin
				for (n, feat) in enumerate(data_int.features)
					lines!(axte, -signal_2[n, :][:]./feat_norm_dict[feat], 1:size(signal_2, 2), color=feat_color_dict[feat],
						linestyle=feat_style_dict[feat])
				end
				axte.topspinevisible=false
				axte.leftspinevisible=false
				axte.rightspinevisible=false
				axte.bottomspinevisible=false
				axte.xgridvisible=false
				axte.ygridvisible=false
				axte.yticklabelsvisible=false
				axte.xticklabelsvisible=false
				axte.yticksvisible=false
				axte.xticksvisible=false
				axte.xlabel = "test (#$test)"
				axte.xlabelsize = 17
			end

			axtr = Axis(gl5[12:14, 3:14])
			begin
				for (n, feat) in enumerate(data_int.features)
					lines!(axtr, 1:size(signal_1, 2), signal_1[n, :]./feat_norm_dict[feat], color=feat_color_dict[feat],
						linestyle=feat_style_dict[feat])
				end
				axtr.topspinevisible=false
				axtr.leftspinevisible=false
				axtr.rightspinevisible=false
				axtr.bottomspinevisible=false
				axtr.xgridvisible=false
				axtr.ygridvisible=false
				axtr.yticklabelsvisible=false
				axtr.xticklabelsvisible=false
				axtr.yticksvisible=false
				axtr.xticksvisible=false
				axtr.ylabel = "train (#$train)"
				axtr.ylabelsize = 17
			end
		end
	end
	Label(f[3, :], text = "", fontsize = 20)
	# save("/Users/joe/Project/PhD/EuroFusion/EUROfusion_ML_2024/Presentations/Multidimensional_qunatity_classification_20_01_25/HYB_BASE_classification.png", F)
	# display(GLMakie.Screen(), F)
	F
end
# plot the final results
let 
	begin
		k = [10, 10]
		train_tss = vcat(Bchoices, Hchoices)
		n = length(train_tss)

		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		shot_dict = Dict([a => n for (n, a) in enumerate(data_int.tok_shots)])
		labelled_ind = [shot_dict[k] for k in labelled_ts]
		train_ind = training_partion(labelled_data, unique(labelled_y), k, S=S)

		X = exp.(-(Array(data_int.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
			exp.(-(Array(data_int.flat_top_cost[!, 2:end]) ./ ft_med).^2)

		K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
		model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed)

		col_dict = Dict("baseline"=>1, "hybrid"=>2, "training"=>3)
		CO_dict = Dict("CO"=> :cross, "C-EH"=> :ltriangle, "C-LH"=> :rtriangle)

		dat = Array(data_int.cosine_cost[labelled_ind, 2:end])
		mat_c = exp.(-(dat[:, Not(labelled_ind)] ./ cos_med).^2)

		dat = Array(data_int.flat_top_cost[labelled_ind, 2:end])
		mat_ft = exp.(-(dat[:, Not(labelled_ind)] ./ ft_med).^2)

		mat = mat_c .* mat_ft
	end

	CairoMakie.activate!()
	f = Figure(size=(1200, 600))
	# gl1 = f[1:2, 1] = GridLayout()
	# let
	# 	axleft = Axis(gl1[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
	# 		limits=(nothing, (0, 1.)))
	# 	axmid = Axis(gl1[1, 2], xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
	# 		limits=(nothing, (0, 1.)))
	# 	axright = Axis(gl1[1, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)))

	# 	ℓ = length(labelled_ind)

	# 	CD = vcat([mat_c[i, i+1:end] for i in 1:ℓ-1]...)
	# 	hist!(axleft, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
	# 	# if σ
	# 	# 	hlines!(axleft, [mat_c[trind, teind]], color=:red)
	# 	# end

	# 	CD = vcat([mat_ft[i, i+1:end] for i in 1:ℓ-1]...)
	# 	hist!(axmid, CD, normalization=:probability, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
	# 	# if σ
	# 	# 	hlines!(axmid, [mat_ft[trind, teind]], color=:red)
	# 	# end
		
	# 	mat = mat_c .* mat_ft
	# 	CD = vcat([mat[i, i+1:end] for i in 1:ℓ-1]...)
	# 	hist!(axright, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
	# 	# if σ
	# 	# 	hlines!(axright, [mat[trind, teind]], color=:red)
	# 	# end
	# end

	gl2 = f[1, 1] = GridLayout()
	begin
		pred_shots = [j for (i, j) in res.ts]
		shot_dict = Dict([a => n for (n, a) in enumerate(data_int.tok_shots)])

		CH = [current_heating("aug", shot) for shot in pred_shots]
		
		training_D = DataFrame(:ts => vcat(Bchoices, Hchoices),
							:predict => vcat(["training" for _ in 1:length(Bchoices)],
								["training" for _ in 1:length(Hchoices)]),
							:conf => vcat([(ts[2], (0.6*cos(n+1))) for (n, ts) in enumerate(Bchoices)],
										[(ts[2], (0.4*cos(n+1))) for (n, ts) in enumerate(Hchoices)]))
		training_D.CH = [current_heating(ts...) for ts in training_D.ts]
		Bn, Hn = length(Bchoices), length(Hchoices)
		tr_base, tr_hyb = training_D[1:Bn, :], training_D[Bn+1:Bn+Hn, :] 
		res.conf = [Tuple.(zip(ts[2], confidence[1, n])) for (n, ts) in enumerate(res.ts)]
		
		res.CH .= CH
		ỹ = res.predict
		
		HYBRID = (@subset res @byrow :predict == "hybrid")
		base = (@subset res @byrow :predict == "baseline")

		ax_base_tr = Axis(gl2[1:3, 1:14])
		ax_pred = Axis(gl2[4:16, 1:14])
		ax_hyb_tr = Axis(gl2[17:19, 1:14])
		
		scatter!(ax_base_tr, training_D.conf[1:length(Bchoices)]; colormap=:tab10, 
			color=3, 
			colorrange = (1, 3),
			strokewidth = 0.1, 
			markersize = 14,
			marker = [CO_dict[label] for label in training_D.CH[1:length(Bchoices)]],
			label = "training"
		)
		ax_base_tr.xticklabelsvisible = false
		ax_base_tr.yticklabelsvisible = false
		hidespines!(ax_base_tr)
		ax_base_tr.limits = (nothing, (-1, 1))
		hidedecorations!(ax_base_tr)

		scatter!(ax_pred, res.conf; colormap=:tab10, 
			color=[col_dict[el] for el in res.predict], 
			colorrange = (1, 3),
			strokewidth = 0.1, 
			markersize = 14,
			marker = [CO_dict[label] for label in res.CH],
			label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ]
		)
		ax_pred.xticks=(10000:10000:40000, ["10,000", "20,000", "30,000", "40,000"])
		
		scatter!(ax_hyb_tr, training_D.conf[length(Bchoices)+1:end]; colormap=:tab10, 
			color=3, 
			colorrange = (1, 3),
			strokewidth = 0.1, 
			markersize = 14,
			marker = [CO_dict[label] for label in training_D.CH[length(Bchoices)+1:end]],
			label = "training"
		)
		ax_hyb_tr.yticklabelsvisible = false
		hidespines!(ax_hyb_tr)
		ax_hyb_tr.limits = (nothing, (-1, 1))
		hidedecorations!(ax_hyb_tr)
		
		linkxaxes!(ax_base_tr, ax_pred, ax_hyb_tr)
		
		pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
		for (n, (ts, coord)) in enumerate(zip(HYBRID.ts, HYBRID.conf))
			text!(ax_pred, coord, text="$(ts[2])", align=pos[mod(n+3,4)+1], fontsize=13, color=:gray50)
		end
		pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
		for (n, (ts, coord)) in enumerate(zip(base.ts, base.conf))
			text!(ax_pred, coord, text="$(ts[2])", align=pos[mod(n+3,4)+1], fontsize=13, color=:gray50)
		end
		pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
		for (n, (ts, coord)) in enumerate(zip(tr_base.ts, tr_base.conf))
			text!(ax_base_tr, coord..., text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13, color=:gray50)
		end
		for (n, (ts, coord)) in enumerate(zip(tr_hyb.ts, tr_hyb.conf))
			text!(ax_hyb_tr, coord..., text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13, color=:gray50)
		end
		
		begin
			group_marker = [MarkerElement(marker = marker, color = :black,
				strokecolor = :transparent, markersize = 12) for marker in [:ltriangle, :rtriangle, :cross]]

			colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
			group_color = [MarkerElement(marker = :circle, color = colors_3[i],
				strokecolor = :transparent, markersize = 12) for i in 1:3]

			Legend(gl2[13:19, 1], [group_marker, group_color],
				[["EH", "LH", "CO"], ["baseline", "hybrid", "training"]],
				["", ""], tellheight=true, framevisible=false, tellwidth=true,
				rowgap=2)
		end
		println("hybrid = ", (@subset res @byrow :predict == "hybrid").CH |> countmap)
		println("baseline = ", (@subset res @byrow :predict == "baseline").CH |> countmap)
	end
	display(GLMakie.Screen(), f)	
end
# compare shots in all 8 parameters
 let
	# 36177, 29772, 13622, 34610, 17282, 15714, 19314, 32135, 34441, 32464, 33407
	x1_lim = 10
	x2_lim = 10
	fig1 = Figure(size=(1000, 400))
	ax1 = Axis(fig1[1:3, 2:4], title="#$train", titlesize=15)
	ax2 = Axis(fig1[4:6, 2:4])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", train)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=:tab10, colorrange=(1, 8), color=n, label=label)
	end
	ax1.limits = ((0., x1_lim), (-0.1, 2.5))
	ax2.limits = ((0., x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	linkxaxes!(ax1, ax2)
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)
	
	ax3 = Axis(fig1[1:3, 5:7], title="#$test", titlesize=15)
	ax4 = Axis(fig1[4:6, 5:7])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", test)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=:tab10, colorrange=(1, 8), color=n, label=label)
	end
	ax3.limits = ((0., x2_lim), (-0.1, 2.5))
	ax4.limits = ((0., x2_lim), (-0.1, 2.5))
	ax3.yticks = 0:0.5:2.5
	ax4.yticks = 0:0.5:2.5
	linkxaxes!(ax3, ax4)
	# fig1[1:2, 8] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	# fig1[4:5, 8] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/shot_comparison.png", fig1)
	display(GLMakie.Screen(), fig1)
end

begin
	present = DataFrame()
	for n in [13, 14]
		dict = Dict(["$i" for i in AUG_shots] .=> -1.)
		
		Bchoices = [("aug", shot) for shot in chosen_baseline_training]	
		para = D[n, :parameters]
		cos_med, ft_med, C = complete_BH[para].hyperparameters

		data_int = data_dtw[para]
		res, confidence, kernel_matrix, model = let
			labelled_data = let
				labelled_shots = vcat(Bchoices, Hchoices)
				labels = vcat(repeat(["baseline"], length(Bchoices)),
					repeat(["hybrid"], length(Hchoices))) 
				d = OrderedDict([ts for ts in labelled_shots] .=> labels)
			end
			classify!(data_int, labelled_data, (cos_med, ft_med, C); interesting="hybrid")
		end
		for (mach, shot) in Hchoices
			dict["$shot"] = 100
		end
		res.conf .= confidence[1, :]
		hyb = (@subset res @byrow :predict == "hybrid")
		for ((mach, shot), conf) in zip(hyb.ts, hyb.conf)
			dict["$shot"] = abs(conf)
		end
		int = DataFrame(dict)
		insertcols!(int, 1, :parameters => para)
		append!(present, int)
	end
end
useful_inds = []
for (n, col) in enumerate(eachcol(present[!, 2:end]))
	if !in(+(col...), [-2, 200])
		push!(useful_inds, n+1)
	end
end
present = present[!, vcat(1, useful_inds)]
present |> vscodedisplay

function increase_decrease(before::DataFrameRow, after::DataFrameRow)
	before, after = Array(before), Array(after)
	[(a_val == -1) ? :xcross : 
		(b_val < a_val ? '↑' : '↓') for (a_val, b_val) in zip(after, before)]
end
function new_point(before::DataFrameRow, after::DataFrameRow)
	before, after = Array(before), Array(after)
	[(b_val == -1 && a_val > 0) ? '😄' : :circle for (a_val, b_val) in zip(after, before)]
end

let
	GLMakie.activate!()
	fig = Figure(size=(1600, 750))
	gl1 = fig[1, 1:3] = GridLayout()

	shots = names(present)[2:end]
	first_half = div(length(shots), 2)
	second_half = length(shots) - div(length(shots), 2)

	colors = vcat(:black, ColorSchemes.colorschemes[:Set2_3]...)
	arrows = [increase_decrease(present[i-1, 2:end], present[i, 2:end]) for i in 2:2][1]
	regular_marker = [new_point(present[i-1, 2:end], present[i, 2:end]) for i in 2:2][1]

	for K in 1:2
		if K == 1
			SHOTS = shots[1:first_half]
			ARROWS = arrows[1:first_half]
			REG = regular_marker[1:first_half] 
			ax = Axis(gl1[K, 1], limits=((-1, first_half+1), (-0.1, 2.2)), 
				xticks=(1:first_half, "#".*SHOTS), ygridvisible=false, xgridvisible=false, 
				xticklabelrotation=π/2, xticklabelsize=13)
		else
			SHOTS = shots[first_half+1:end]
			ARROWS = arrows[first_half+1:end]
			REG = regular_marker[first_half+1:end]
			ax = Axis(gl1[K, 1], limits=((-1, second_half+1), (-0.1, 2.2)), 
				xticks=(1:second_half, "#".*SHOTS), ygridvisible=false, xgridvisible=false, 
				xticklabelrotation=π/2, xticklabelsize=13)
		end

		begin
			for (n, shot) in enumerate(SHOTS)
				col = present[!, shot]
				for (k, val) in enumerate(col)
					# if k <= K
						if in(val, [-1, 100])
							continue
						elseif k == 1
							scatter!(ax, n, val, color=colors[k], marker=ARROWS[n], markersize=14, strokewidth=0.1)	
						elseif (k == 2)
							scatter!(ax, n, val, color=colors[k], marker=REG[n])
						# else
						# 	scatter!(ax, n, val, color=(colors[k], transp[k]), label=present.parameters[k])
						end
					# else 
					# 	continue
					# end

				end
			end

			begin
				text!((0, 1.9), text=L"\beta_p", align=(:left, :baseline), fontsize=20, color=(:black))
				text!((0, 1.7), text=L"\beta_p \,\,\, n_{GW}", align=(:left, :baseline), fontsize=20, color=(colors[2]), strokewidth=0.3, strokecolor=colors[2])
			end
		end
		if K == 1
			group_marker = [MarkerElement(marker = marker, color = color,
				strokecolor = :transparent,
				markersize = 12) for (marker, color) in zip([:circle, '↑', '↓', :xcross, '😄'],
					[colors[2], :black, :black, :black, colors[2]])]

			Legend(fig[1, 3],
				group_marker,
				["hybrid prediction", "increased confidence", "decreased confidence", "removed since model 1", "added by model 2"],
				framevisible=false, labelsize=20, halign=:right, valign=:top)
		end
	end
	linkxaxes!(contents(gl1))
	# save("/Users/joe/Project/PhD/EuroFusion/EUROfusion_ML_2024/Presentations/Multidimensional_qunatity_classification_20_01_25/difference_in_βp_nGW.png", fig)
	fig
end














































let
	GLMakie.activate!()
	fig = Figure(size=(1600, 750))
	gl1 = fig[1, 1:3] = GridLayout()
	# gl2 = fig[1, 4] = GridLayout()
	shots = names(present)[2:end]
	for (K, N) in zip(1:2, [13, 14])
		if K == 4
			ax = Axis(gl1[K, 1], limits=((-4, length(shots)+1), (-0.1, 1.8)), 
				xticks=(1:length(shots), "#".*shots), ygridvisible=false, xgridvisible=false, 
				xticklabelrotation=π/2, xticklabelsize=13)
		else
			ax = Axis(gl1[K, 1], limits=((-4, length(shots)+1), (-0.1, 1.8)), ygridvisible=false)
			hidexdecorations!(ax)
		end
		transp = reverse([(1/k)^2 for k in 1:K])
		begin
			colors = vcat(:black, ColorSchemes.colorschemes[:Set2_3]...)
			arrows = [increase_decrease(present[i-1, 2:end], present[i, 2:end]) for i in 2:2]
			regular_marker = [new_point(present[i-1, 2:end], present[i, 2:end]) for i in 2:2]

			for (n, shot) in enumerate(shots)
				col = present[!, shot]
				for (k, val) in enumerate(col) 
					if k <= K
						if in(val, [-1, 100])
							continue
						elseif k == K-1
							scatter!(ax, n, val, color=colors[k], label=present.parameters[k], marker=arrows[k][n], markersize=14, strokewidth=0.1)	
						elseif (k == K && K != 1)
							scatter!(ax, n, val, color=(colors[k], transp[k]), label=present.parameters[k], marker=regular_marker[k-1][n])
						else
							scatter!(ax, n, val, color=(colors[k], transp[k]), label=present.parameters[k])
						end
					else 
						continue
					end
				end
			end

			begin
				text!((-3.5, 1.55), text=L"\beta_p", align=(:left, :baseline), fontsize=17, color=(:black, 1))
				text!((-3.5, 1.32), text=L"I_p", align=(:left, :baseline), fontsize=17, color=(:black, 1))
				text!((-3.5, 1.09), text=L"\ell_i", align=(:left, :baseline), fontsize=17, color=(:black, (K > 1 ? 0. : 1)))
				text!((-3.5, 0.86), text=L"n_{GW}", align=(:left, :baseline), fontsize=17, color=(:black, (K > 2 ? 0. : 1)))
				text!((-3.5, 0.63), text=L"P_{ECRH}", align=(:left, :baseline), fontsize=17, color=(:black, 1))
				text!((-3.5, 0.41), text=L"P_{ICRH}", align=(:left, :baseline), fontsize=17, color=(:black, (K > 3 ? 0. : 1)))
				text!((-3.5, 0.18), text=L"P_{NBI}", align=(:left, :baseline), fontsize=17, color=(:black, 1))
				text!((-3.5, -0.05), text=L"Q_{95}", align=(:left, :baseline), fontsize=17, color=(:black, 1))
			end
		end
		begin
			para = D[n, :parameters]
			RES = complete_BH[para]
			cos_med, ft_med, C = RES.hyperparameters

			gl2K = gl2[K, 1] = GridLayout()
			axleft = Axis(gl2K[1, 1], xticks=(1:0.2:1, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
				limits=(nothing, (0, 1.)), title="σ=$(round(cos_med, digits=1))")
			axmid = Axis(gl2K[1, 2], xticks=(1:0.2:1, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
				limits=(nothing, (0, 1.)), title="σ=$(round(ft_med, digits=1))")
			axright = Axis(gl2K[1, 3], xticks=(1:1:1, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)))
				transp = reverse([(1/k)^2 for k in 1:K])

			linkyaxes!(axleft, axright, axmid)
		end
		# begin
		# 	ℓ = 100
		# 	plot_D = (@subset RES.results @byrow :parameters == para)

		# 	boxplot!(axleft, ones(ℓ), plot_D.cos_hist_lb, color=:gray70)
		# 	boxplot!(axleft, ones(ℓ), plot_D.cos_hist_ub, color=:gray50)
		# 	boxplot!(axleft, ones(ℓ), plot_D.cos_mode, color=:red)

		# 	boxplot!(axmid, ones(ℓ), plot_D.ft_hist_lb, color=:gray70)
		# 	boxplot!(axmid, ones(ℓ), plot_D.ft_hist_ub, color=:gray50)
		# 	boxplot!(axmid, ones(ℓ), plot_D.ft_mode, color=:red)

		# 	boxplot!(axright, ones(ℓ), plot_D.kernel_hist_lb, color=:gray70)
		# 	boxplot!(axright, ones(ℓ), plot_D.kernel_hist_ub, color=:gray50)
		# 	boxplot!(axright, ones(ℓ), plot_D.kernel_mode, color=:red)
		# end
	end
	linkxaxes!(contents(gl1))
	# save("/Users/joe/Project/PhD/EuroFusion/EUROfusion_ML_2024/Presentations/Multidimensional_qunatity_classification_20_01_25/choosing_parameters.png", fig)
	fig
end

wine_tasting = [15632, 16933, 18439, 18884, 19112, 28882, 8042]






















# data_CEL = let
#     CEL_features = ["PNBI"]
#     CEL_tss = tok_shots((@subset which(CEL_features) @byrow in(:tok, ["aug"])))
#     DTW_hyp_1(CEL_tss, CEL_features, 10, L=100)
# end
p2 = vcat([["IP", "PNBI"]], vcat.("PNBI", ["PECRH", "PICRH", "BETAPOL", "Q95", "NGW", "LI"]))
for para in ["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "NGW", "LI"]
	data_CEL = data_bhind[*([para], delim="_")]
	CLn, CEn, COn = 20, 20, 10
	for S in 1:3
		CLchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false))
		CEchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false))
		COchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false))

		k = [6,6,6]
		
		res, confidence, kernel_matrix, (cos_med, ft_med), ACC = let
			labelled_data = let
				labelled_shots = vcat(CLchoices, CEchoices, COchoices)
				labels = vcat(repeat(["C-LH"], length(CLchoices)),
					repeat(["C-EH"], length(CEchoices)),
					repeat(["CO"], length(COchoices))) 
				d = OrderedDict([ts for ts in labelled_shots] .=> labels)
			end
			try
				levels = ["C-LH", "C-EH", "CO"]
				hyper_parameter_search(data_CEL, labelled_data, k, interesting="CO", N=60, metric="BalancedAccuracy()")
			catch
				continue
			end
		end

		train_tss = vcat(CLchoices, CEchoices, COchoices)
		n = length(train_tss)
		ind = [data_CEL.shot_dict[ts] for ts in train_tss]

		dat = data_CEL.cosine_cost
		mat = Array(dat[ind, ind.+1])
		CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

		cos_mode = mode(CD)
		cos_ext = extrema(CD)

		dat = data_CEL.flat_top_cost
		mat = Array(dat[ind, ind.+1])
		CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

		ft_mode = mode(CD)
		ft_ext = extrema(CD)

		shots = [j for (i, j) in res.ts]
		labels = (@subset all_2D_features @byrow begin
			:tok == "aug"
			in(:shot, shots)
		end).current_heating

		acc, b_acc = Accuracy()(res.predict, labels), BalancedAccuracy()(res.predict, labels)

		append!(meta_D, 
			DataFrame(:parameters => *([para], delim="_"), 
				:random_seed => S,
				:metric => ACC,
				:acc_D => acc,
				:Bacc_D => b_acc,
				:COS => cos_med,
				:cos_hist_mode => cos_mode,
				:cos_hist_ext => cos_ext,
				:FT => ft_med,
				:ft_hist_mode => ft_mode,
				:ft_hist_ext => ft_ext
				))
	end
end

meta_D |> vscodedisplay

Random.seed!(123)
CLn = 20
CEn = 20
COn = 10
CLchoices = sort(sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false))
CEchoices = sort(sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false))
COchoices = sort(sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false))

using CairoMakie
let
	CairoMakie.activate!()
	see = [rand(COchoices), rand(CEchoices), rand(CLchoices)]
	println(see)
	f = Figure()
	a0 = Axis(f[1, 1])
	a1 = Axis(f[1, 2])

	hidedecorations!(a0)
	hidedecorations!(a1)

	
	for (n, (el, label)) in enumerate(zip(see, ["overshoot", "early", "late"]))
		see_pd = data_CEL.profile_data[el]
		lines!(a0, see_pd[2, :], see_pd[1, :]; colormap=:tab10,
			colorrange=(1, 3),
			color=n,
			label=label
		)
	end
	for (n, (el, label)) in enumerate(zip(see, ["overshoot", "early", "late"]))
		see_pd = data_CEL.flat_top_data[el]
		lines!(a1, see_pd[2, :], see_pd[1, :]; colormap=:tab10,
			colorrange=(1, 3),
			color=n,
		)
	end
	axislegend(a0, framevisible=false)
	f
end

data_CEL = data_ind["PNBI"]
res, confidence, kernel_matrix, (cos_med, ft_med)= let
	labelled_data = let
	    labelled_shots = vcat(CLchoices, CEchoices, COchoices)
	    labels = vcat(repeat(["C-LH"], length(CLchoices)),
	        repeat(["C-EH"], length(CEchoices)),
	        repeat(["CO"], length(COchoices))) 
	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	hyper_parameter_search(data_CEL, labelled_data, [6,6,6], interesting="CO", N=60, metric="BalancedAccuracy()")
end
res, confidence, kernel_matrix, model= let
	labelled_data = let
	    labelled_shots = vcat(CLchoices, CEchoices, COchoices)
	    labels = vcat(repeat(["C-LH"], length(CLchoices)),
			repeat(["C-EH"], length(CEchoices)),
			repeat(["CO"], length(COchoices))) 
	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	classify!(data_CEL, labelled_data, (14.47, 1262.46), interesting="CO")
end

S = 123
cos_med = 11.711700000000002
ft_med = 1407.3832499999999
let
	colors = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
	color_dict = Dict("C-LH" => colors[1], "C-EH" => colors[2], "CO" => colors[3])
	k = [6,6,6]

	train_tss = vcat(CLchoices, CEchoices, COchoices)
	n = length(train_tss)

	labelled_data = let
		labelled_shots = vcat(CLchoices, CEchoices, COchoices)
		labels = vcat(repeat(["C-LH"], CLn),
			repeat(["C-EH"], CEn),
			repeat(["CO"], COn)) 
		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end	
	labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data_CEL.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

	train_ind = training_partion(labelled_data, unique(labelled_y), k, S=24)

	dat = Array(data_CEL.cosine_cost[labelled_ind[train_ind], 2:end])
	mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ cos_med).^2)

	dat = Array(data_CEL.flat_top_cost[labelled_ind[train_ind], 2:end])
	mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ ft_med).^2)

	mat = mat_c .* mat_ft

	f = Figure();
	a1 = Axis(f[1, 1])
	# a1 = Axis(f[1:3, 1:20], title="", ylabel="", ylabelsize=20, yticks=1:50:500, ylabelpadding=20)
	# a2 = Axis(f[1:3, 21:25], yticklabelsvisible=false, yticksvisible=false, title="")
	he = heatmap!(a1, mat, colorscale=exp10, colormap=:oslo)

	cntmap = countmap(labelled_y[Not(train_ind)])
	lines!(a1, [k[1]+0.5], [0.5, cntmap["C-LH"]+cntmap["C-EH"]+0.5], color=:white, linewidth=4)
	lines!(a1, [k[1]+k[2]+0.5], [cntmap["C-LH"], length(labelled_y[Not(train_ind)])].+0.5, color=:white, linewidth=4)

	lines!([0, k[1]+k[2]].+0.5, [cntmap["C-LH"]+0.5], color=:white, linewidth=4)
	lines!([k[1], length(labelled_y[train_ind])].+0.5, [cntmap["C-LH"]+cntmap["C-EH"]+0.5], color=:white, linewidth=4)

	# lines!([k[1], k[1]].+0.5, [0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)
	# lines!([k[1]+k[2], k[1]+k[2]].+0.5, [k[1]+0.5, +(k...)+0.5], color=:white, linewidth=4)
	# lines!([0.5, k[1]+k[2]+0.5], [k[1], k[1]].+0.5, color=:white, linewidth=4)
	# lines!([k[1]+0.5, +(k...)+0.5], [k[1]+k[2]+0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)

	Colorbar(f[1, 2], he, ticks=0:0.1:1)
	# a1.xticks = (1:length(Bchoices), ["#$j" for (i,j) in Bchoices])
	a1.xticks = (1:length(labelled_ts[train_ind]), ["#$j" for (i, j) in labelled_ts[train_ind]])
	a1.xticklabelrotation = π/5
	a1.yticks = (1:length(labelled_ts[Not(train_ind)]), ["#$j" for (i, j) in labelled_ts[Not(train_ind)]])
	
	# a1.xticklabelrotation = π/2
	# a2.xticklabelrotation = π/2
	
	# a1.xticklabelsize = 15
	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/BH_kernel_matrix.png", f)
	f
end
let
	# CairoMakie.activate!()
	k = [CLn, CEn, COn]
	CLchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false))
	CEchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false))
	COchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false))
	
	train_tss = vcat(CLchoices, CEchoices, COchoices)
	n = length(train_tss)

	labelled_data = let
		labelled_shots = vcat(CLchoices, CEchoices, COchoices)
		labels = vcat(repeat(["C-LH"], length(CLchoices)),
			repeat(["C-EH"], length(CEchoices)),
			repeat(["CO"], length(COchoices))) 
		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data_CEL.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

	train_ind = training_partion(labelled_data, unique(labelled_y), floor.(Int, 0.6.*k))
	n = sum(train_ind)

	fig = Figure();
	ax = Axis(fig[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20)
	ax1 = Axis(fig[1, 2], xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20)
	ax2 = Axis(fig[1, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20)

	dat = Array(data_CEL.cosine_cost[labelled_ind[train_ind], 2:end])
	mat_c = exp.(-(dat[:, labelled_ind[train_ind]] ./ cos_med).^2)
	CD = vcat([mat_c[i, i+1:end] for i in 1:n-1]...)
	hist!(ax, CD, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
	# hlines!(ax, [COS_med], xmin=0.05, xmax=0.95, color=:red)

	dat = Array(data_CEL.flat_top_cost[labelled_ind[train_ind], 2:end])
	mat_ft = exp.(-(dat[:, labelled_ind[train_ind]] ./ ft_med).^2)
	CD = vcat([mat_ft[i, i+1:n] for i in 1:n-1]...)
	hist!(ax1, CD, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
	# # hlines!(ax1, [FT_med], xmin=0.05, xmax=0.95, color=:red)
	
	mat = mat_c .* mat_ft
	CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
	hist!(ax2, CD, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
	# # # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/cost_spread_BH.png", fig)
	fig
end
let
	# GLMakie.activate!()
	shots = [j for (i, j) in res.ts]
	labels = (@subset all_2D_features @byrow begin
	    :tok == "aug"
	    in(:shot, shots)
	end).current_heating
	res.label .= labels
	res.conf = Tuple.(zip(confidence[1, :], confidence[2, :], confidence[3, :]))
	curr_over = (@subset res @byrow :predict == "CO")

	acc, b_acc = Accuracy()(res.predict, labels), BalancedAccuracy()(res.predict, labels)
	println("accuracy = $(acc), balanced accuracy = $(b_acc)")

	col_dict = Dict("CO"=>1, "C-EH"=>2, "C-LH"=>3)
    ỹ = res.predict
    y = res.label

	f = Figure()
	a = Axis3(f[1,1])
   
    scatter!(confidence[1, :], confidence[2, :], confidence[3, :]; colormap=:tab10, 
        color=[col_dict[el] for el in ỹ],
		strokewidth = 0.1,
        label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
		markersize = [(i == j ? 20 : 10) for (i, j) in zip(y, ỹ)],
        marker = [(i == j ? '∘' : (:utriangle)) for (i, j) in zip(y, ỹ)]
    )
	
    pos = [(:right, :bottom), (:left, :bottom), (:right, :top), (:left, :top)]
    Random.seed!(123)
    shuffle!(curr_over)
    for (n, (ts, coord)) in enumerate(zip(curr_over.ts, curr_over.conf))
        text!(a, coord..., text="$(ts[2])", align=pos[mod(n,4)+1], fontsize=7)
    end

	# extra = (@subset res @byrow in(:ts, [("aug", shot) for shot in [27930]]))
	# for (n, (ts, coord)) in enumerate(zip(extra.ts, extra.conf))
    #     text!(a, coord, text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=7)
    # end

    axislegend(a, unique=true, merge=true, position=:lb)
    a.xgridvisible = false
    a.ygridvisible = false
	a.zgridvisible = false
    
	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/CO_EH_LH_classification.png", f)
	f
end

let
	train_tss = vcat(CLchoices, CEchoices, COchoices)
	n = length(train_tss)
	ind = [data_CEL.shot_dict[ts] for ts in train_tss]

	fig = Figure();
	ax = Axis(fig[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20)
	ax1 = Axis(fig[1, 2], xticks=(0.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{I_p^{80%}}", L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20)

	dat = data_CEL.cosine_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))
	hist!(ax, CD, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
	println("mode = ", mode(CD), "\n", 
		"extrema = ", extrema(CD))
	hlines!(ax, [cos_med], xmin=0.05, xmax=0.95, color=:red)

	dat = data_CEL.flat_top_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))
	hist!(ax1, CD, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
	println("mode = ", mode(CD), "\n", 
		"extrema = ", extrema(CD))
	hlines!(ax1, [ft_med], xmin=0.05, xmax=0.95, color=:red)
	
	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/cost_spread.png", fig)
	fig
end

let
    shot = 11222
    shot_dict = Dict(res.ts .=> 1:length(res.ts))
	closeness = kernel_matrix[:, shot_dict[("aug", shot)]]
	D = DataFrame(:label => vcat(["C-LH" for _ in 1:length(CLchoices)], ["C-EH" for _ in 1:length(CEchoices)], ["CO" for _ in 1:length(COchoices)]), 
                    :ts => vcat(CLchoices, CEchoices, COchoices), 
                    :dist => closeness)
	
	(@subset res @byrow begin
		:predict == "C-LH"
		:label == "C-EH"
	end), sort(D, :dist, rev=true)
end

let
	train = 23074
	test = 19112
	# 36177, 29772, 13622, 34610, 17282, 15714, 19314, 32135, 34441, 32464, 33407, 34841 

	tr_ind = findall(i -> in(i, [("aug", train)]), data_CEL.tok_shots)[1]
	te_ind = findall(i -> in(i, [("aug", test)]), data_CEL.tok_shots)[1]

	CC = round(data_CEL.cosine_cost[tr_ind, te_ind+1], digits=2)
	FTC = round(data_CEL.flat_top_cost[tr_ind, te_ind+1], digits=2)
	
	CCexp = round(exp(-CC^2 / cos_med^2), digits=2)
	FTCexp = round(exp(-FTC^2 / ft_med^2), digits=2)
	
	shot = train
	x1_lim = 6.8
	x2_lim = 6.5
	fig1 = Figure(size=(1900, 900))
	ax1 = Axis(fig1[1:3, 2:4], title="#$shot", titlesize=25)
	ax2 = Axis(fig1[4:6, 2:4])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=2.5, colormap=:dracula, colorrange=(1, 8), color=n, label=label)
	end
	ax1.limits = ((0.4, x1_lim), (-0.1, 2.5))
	ax2.limits = ((0.4, x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	text!(ax1, (0.5, 2.2), text = "cosine cost = $(CCexp)", fontsize=25)
	text!(ax1, (0.5, 1.6), text = "flat top cost = $(FTCexp)", fontsize=25)
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)
	
	shot = test
	ax3 = Axis(fig1[1:3, 5:7], title="#$shot", titlesize=25)
	ax4 = Axis(fig1[4:6, 5:7])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=2.5, colormap=:dracula, colorrange=(1, 8), color=n, label=label)
	end
	ax3.limits = ((0.4, x2_lim), (-0.1, 2.5))
	ax4.limits = ((0.4, x2_lim), (-0.1, 2.5))
	ax3.yticks = 0:0.5:2.5
	ax4.yticks = 0:0.5:2.5
	# fig1[1:2, 8] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	# fig1[4:5, 8] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/shot_comparison.png", fig1)
	fig1
end





let 
	f = Figure(size=(1000, 1200));
	gl1 = GridLayout()
	gl1[1, 1] = Axis(f)
	let
		k = [10, 10]

		train_tss = vcat(Bchoices, Hchoices)
		n = length(train_tss)

		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])
		labelled_ind = [shot_dict[k] for k in labelled_ts]

		train_ind = training_partion(labelled_data, unique(labelled_y), k, S=10)

		dat = Array(data_BH.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ COS_med).^2)

		dat = Array(data_BH.flat_top_cost[labelled_ind[train_ind], 2:end])
		mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ FT_med).^2)

		mat = mat_c .* mat_ft
		
		# a1 = Axis(f[1:3, 1:20], title="", ylabel="", ylabelsize=20, yticks=1:50:500, ylabelpadding=20)
		# a2 = Axis(f[1:3, 21:25], yticklabelsvisible=false, yticksvisible=false, title="")
		he = heatmap!(contents(gl1)[1], mat, colorscale=exp10, colormap=:oslo)

		cntmap = countmap(labelled_y[Not(train_ind)])
		lines!(contents(gl1)[1], [k[1]+0.5], [0.5, cntmap["baseline"]+cntmap["hybrid"]+0.5], color=:white, linewidth=4)
		# lines!(a1, [k[1]+k[2]+0.5], [cntmap["C-LH"], length(labelled_y[Not(train_ind)])].+0.5, color=:white, linewidth=4)

		lines!(contents(gl1)[1], [0, k[1]+k[2]].+0.5, [cntmap["baseline"]+0.5], color=:white, linewidth=4)
		# lines!([k[1], length(labelled_y[train_ind])].+0.5, [cntmap["C-LH"]+cntmap["C-EH"]+0.5], color=:white, linewidth=4)

		# lines!([k[1], k[1]].+0.5, [0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)
		# lines!([k[1]+k[2], k[1]+k[2]].+0.5, [k[1]+0.5, +(k...)+0.5], color=:white, linewidth=4)
		# lines!([0.5, k[1]+k[2]+0.5], [k[1], k[1]].+0.5, color=:white, linewidth=4)
		# lines!([k[1]+0.5, +(k...)+0.5], [k[1]+k[2]+0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)

		# a1.xticks = (1:length(Bchoices), ["#$j" for (i,j) in Bchoices])
		contents(gl1)[1].yticks = (1:length(labelled_ts[Not(train_ind)]), ["#$j" for (i, j) in labelled_ts[Not(train_ind)]])
		contents(gl1)[1].xticks = (1:length(labelled_ts[train_ind]), ["#$j" for (i, j) in labelled_ts[train_ind]])
		contents(gl1)[1].xticklabelrotation = π/5
		
		# a1.xticklabelrotation = π/2
		# a2.xticklabelrotation = π/2
		
		# a1.xticklabelsize = 15
		# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/BH_kernel_matrix.png", f)
	end
	gl2 = GridLayout()
	gl2[1, 1:3] = [Axis(f, xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20),
		Axis(f, xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20),
		Axis(f, xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20)]
	let
		# CairoMakie.activate!()
		k = [10, 10]
		
		train_tss = vcat(Bchoices, Hchoices)
		n = length(train_tss)

		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])
		labelled_ind = [shot_dict[k] for k in labelled_ts]

		train_ind = training_partion(labelled_data, unique(labelled_y), k)
		n = sum(train_ind)

		# fig = Figure();
		# ax2 = Axis(f[1, 2])
		# ax2 = Axis(fig[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20)
		# ax1 = Axis(fig[1, 2], xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20)
		# ax2 = Axis(fig[1, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20)

		dat = Array(data_BH.cosine_cost[labelled_ind[train_ind], 2:end])
		mat_c = exp.(-(dat[:, labelled_ind[train_ind]] ./ COS_med).^2)
		CD = vcat([mat_c[i, i+1:end] for i in 1:n-1]...)
		hist!(gl2[1, 1], CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
		# hlines!(ax, [COS_med], xmin=0.05, xmax=0.95, color=:red)

		dat = Array(data_BH.flat_top_cost[labelled_ind[train_ind], 2:end])
		mat_ft = exp.(-(dat[:, labelled_ind[train_ind]] ./ FT_med).^2)
		CD = vcat([mat_ft[i, i+1:n] for i in 1:n-1]...)
		hist!(gl2[1, 2], CD, normalization=:probability, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
		# # hlines!(ax1, [FT_med], xmin=0.05, xmax=0.95, color=:red)
		
		mat = mat_c .* mat_ft
		CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
		hist!(gl2[1, 3], CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
		# # # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/cost_spread_BH.png", fig)
		# fig
	end
	gl3 = f[2:3, 1:2] = GridLayout()
	# gl3[1:3, 1] = [Axis(f) for _ in 1:3]
	# rowsize!(gl3, 1, Relative(0.1))
	# rowsize!(gl3, 3, Relative(0.1))
	let
		pred_shots = [j for (i, j) in RES.ts]
		shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])

		CH = [current_heating("aug", shot) for shot in pred_shots]
		
		col_dict = Dict("baseline"=>1, "hybrid"=>2, "training"=>3)
		CO_dict = Dict("CO"=> :cross, "C-EH"=> :ltriangle, "C-LH"=> :rtriangle)

		training_D = DataFrame(:ts => vcat(Bchoices, Hchoices),
							:predict => vcat(["training_baseline" for _ in 1:length(Bchoices)],
								["training_hybrid" for _ in 1:length(Hchoices)]),
							:conf => vcat([(ts[2], (0.6*cos(n+1))) for (n, ts) in enumerate(Bchoices)],
										[(ts[2], (0.55*cos(n+1))) for (n, ts) in enumerate(Hchoices)]))
		training_D.CH = [current_heating(ts...) for ts in training_D.ts]

		RES.conf = [Tuple.(zip(ts[2], CONFIDENCE[1, n])) for (n, ts) in enumerate(RES.ts)]
		
		# RES_updated = vcat(RES[!, [:ts, :predict, :conf]], training_D)
		RES.CH .= CH
		ỹ = RES.predict
		
		HYBRID = (@subset RES @byrow :predict == "hybrid")[1:2:end, :]
		CO = (@subset RES @byrow begin
			:CH == "CO"
			:predict == "baseline"
		end)

		axtop = Axis(gl3[1:3, 1:12])
		axmid = Axis(gl3[4:16, 1:12], ylabel="confidence", xlabel="shot #")
		axbot = Axis(gl3[17:19, 1:12])
		a_sp = scatter!(axmid, RES.conf; colormap=:tab10, 
			color=[col_dict[el] for el in RES.predict], 
			colorrange = (1, 3),
			strokewidth = 0.1, 
			markersize = 18,
			marker = [CO_dict[label] for label in RES.CH],
			label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ])
		axtop.xticklabelsvisible = false
		axtop.yticklabelsvisible = false
		hidespines!(axtop)
		axtop.limits = (nothing, (-0.9, 0.9))
		hidedecorations!(axtop)
		axmid.xticks=(10000:10000:40000, ["10,000", "20,000", "30,000", "40,000"])
		axbot.yticklabelsvisible = false
		hidespines!(axbot)
		axbot.limits = (nothing, (-0.9, 0.9))
		hidedecorations!(axbot)
		linkxaxes!(axtop, axmid, axbot)
		
		a0_sp = scatter!(axtop, training_D.conf[1:length(Bchoices)]; colormap=:tab10, 
			color=3, 
			colorrange = (1, 3),
			strokewidth = 0.1, 
			markersize = 18,
			marker = [CO_dict[label] for label in training_D.CH[1:length(Bchoices)]],
			label = "training")

		scatter!(axbot, training_D.conf[length(Bchoices)+1:end]; colormap=:tab10, 
			color=3, 
			colorrange = (1, 3),
			strokewidth = 0.1, 
			markersize = 18,
			marker = [CO_dict[label] for label in training_D.CH[length(Bchoices)+1:end]],
			label = "training")

		pos = [(:left, :top), (:right, :bottom), (:left, :bottom), (:right, :top)]
		for (n, (ts, coord)) in enumerate(zip(HYBRID.ts, HYBRID.conf))
			text!(axmid, coord, text="$(ts[2])", align=pos[mod(n+0,4)+1], fontsize=13)
		end
		pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
		for (n, (ts, coord)) in enumerate(zip(CO.ts, CO.conf))
			text!(axmid, coord, text="$(ts[2])", align=pos[mod(n+0,4)+1], fontsize=13)
		end
		pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
		tr_D_hy = @subset training_D @byrow :predict == "training_hybrid"
		for (n, (ts, coord)) in enumerate(zip(tr_D_hy.ts, tr_D_hy.conf))
			text!(axbot, coord, text="$(ts[2])", align=pos[mod(n+0,4)+1], fontsize=13)
		end
		pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
		tr_D_CO = @subset training_D @byrow :predict == "training_baseline"
		for (n, (ts, coord)) in enumerate(zip(tr_D_CO.ts, tr_D_CO.conf))
			text!(axtop, coord, text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13)
		end

		group_marker = [MarkerElement(marker = marker, color = :black,
			strokecolor = :transparent,
			markersize = 12) for marker in [:ltriangle, :rtriangle, :cross]]

		Legend(f[2:3, 1:2][14, 1],
			group_marker,
			["EH", "LH", "CO"],
			framevisible=false, labelsize=20, halign=:left, valign=:top)

		colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
		group_marker = [MarkerElement(marker = :circ, color = colors_3[i],
			strokecolor = :transparent,
			markersize = 12) for i in 1:3]

		Legend(f[2:3, 1:2][17, 1],
			group_marker,
			["baseline", "hybrid", "training"],
			framevisible=false, labelsize=20, halign=:left, valign=:top)

		# # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/Baseline_Hybrid_classified.png", f)
		
		# # println("hybrid = ", (@subset RES @byrow :predict == "hybrid").CH |> countmap)
		# # println("baseline = ", (@subset RES @byrow :predict == "baseline").CH |> countmap)
		
		# f
	end
	
	f.layout[1, 1] = gl1
	f.layout[1, 2] = gl2
	# f.layout[2:3, 1:2] = gl3
	Label(f[0, :], text = "", fontsize = 20)
	
	f
end










# all_IB = [7643, 8127, 29619, 29620, 29621, 29622, 29623, 29636, 29953, 29956, 29957, 29958, 29962, 29963, 29964, 29965, 29976, 29977, 29979, 29980, 29981, 32103, 32104, 32105, 32106, 32109, 32110, 32111, 32112, 32113, 32114, 32115, 32120, 32121, 32122, 32123, 32463, 32464, 32471, 32472, 32474, 32475, 32489, 33370, 33371, 33374, 33375, 33376, 33377, 33406, 33407, 33409, 33427, 34441, 34442, 34443, 34444, 34445, 34446, 34449, 34450, 34454, 34489, 34490, 34840, 35975, 36177, 37081, 37094, 37096, 39124, 39125, 39126, 39128, 39129, 40411, 40851, 40410, 37093, 37082, 36178, 37080, 36176, 34842, 40206, 40207, 40202, 40203, 34838, 34839, 34841, 34844
# ]
# 23259, 23971
chosen_baseline_training = [6483, 7558, 7643, 8007, 8093, 10729, 11126, 13603, 15714, 17389, 20280,  
	25719, 25757, 29623, 32225, 33303, 34489, 34490, 35975, 40203, 40411
]
Bchoices = [("aug", shot) for shot in chosen_baseline_training]
Hchoices = [("aug", shot) for shot in [11190, 16688, 16736, 18046, 18869, 18880, 19314, 25764, 26338, 26913, 27930, 34769, 34770, 
	34774, 34775, 34776, 36443, 36408]]

Bn = length(Bchoices)
Hn = length(Hchoices)

# data_bhind = data_ind
# data_bh1 = Dict()
# for (key, val) in data_bhind
# 	new_key = *(sort(split(key, "_")), delim="_")
# 	data_bh1[new_key] = val
# end
# sort(data_bh1 |> keys |> collect)

# = let
#     dict = Dict()
#     for feat in ["IP", "PNBI", "PECRH", "PICRH", "BETAPOL", "Q95", "NGW", "LI"]
#         tss = tok_shots((@subset which(feat) @byrow in(:tok, ["aug"])))
#         dict[feat] = DTW_hyp_1(tss, [feat], 10, L=100)
#     end
#     dict
# end

# p1 = [["IP"], ["PNBI"], ["PICRH"], ["PECRH"], ["BETAPOL"], ["Q95"], ["NGW"], ["LI"]]
# p3 = vcat.([["BETAPOL", "PNBI"]], ["IP", "PECRH", "PICRH", "Q95", "NGW", "LI"])
# p4 = vcat.([["IP", "PNBI", "BETAPOL"]], ["PECRH", "PICRH", "Q95", "NGW", "LI"])
# p5 = vcat.([["IP", "PNBI", "BETAPOL", "NGW"]], ["PECRH", "PICRH", "Q95", "LI"])
# p6 = vcat.([["IP", "PNBI", "PECRH", "BETAPOL", "NGW"]], ["PICRH", "Q95", "LI"])
for feat in [["IP", "PNBI", "PECRH", "PICRH", "Q95", "NGW", "BETAPOL"]]
	if in(*(sort(feat), delim="_"), collect(keys(data_bhind)))
		continue
	end
    tss = tok_shots((@subset which(feat) @byrow in(:tok, ["aug"])))
    data_bhind[*(sort(feat), delim="_")] = DTW_hyp_1(tss, feat, 10, L=100)
end
data_bhind

# meta_D_BH = DataFrame()
# for (n, (para, cos, ft)) in enumerate(eachrow(meta_D_BH[!, [:parameters, :COS, :FT]]))
# 	for S in 1:50	
# 		labelled_data = let
# 			labelled_shots = vcat(Bchoices, Hchoices)
# 			labels = vcat(repeat(["baseline"], length(Bchoices)),
# 				repeat(["hybrid"], length(Hchoices))) 
# 			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
# 		end
# 		labelled_ts = collect(labelled_data |> keys)
# 		labelled_y = collect(labelled_data |> values)
# 		shot_dict = Dict([a => n for (n, a) in enumerate(data_bhind[para].tok_shots)])
# 		labelled_ind = [shot_dict[k] for k in labelled_ts]

# 		shot_dict = Dict([a => n for (n, a) in enumerate(data_bhind[para].tok_shots)])
# 		labels = collect(labelled_data |> values) |> unique

# 		train_ind = training_partion(labelled_data, labels, k, S=S)

# 		X = exp.(-(Array(data_bhind[para].cosine_cost[!, 2:end]) ./ cos).^2) .* 
# 			exp.(-(Array(data_bhind[para].flat_top_cost[!, 2:end]) ./ ft).^2)

# 		K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
# 		model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed)

# 		KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
# 		ỹ, _ = svmpredict(model, KK)
# 		# ỹ, _ = KNN(KK, labelled_y[train_ind])
# 		meta_D_BH[n, :Bacc_D] += BalancedAccuracy()(ỹ, labelled_y[Not(train_ind)])
# 	end
# end
# meta_D_BH.Bacc_D ./= 50
# meta_D_BH.Bacc_D .= 0.
# meta_D_BH |> vscodedisplay

for (para, data) in data_bhind
	if in(para, meta_D_BH.parameters)
		continue
	end
	data_CEL = data
	k = [10, 10]

	res, confidence, kernel_matrix, (cos_med, ft_med), ACC = let
		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		try
			levels = ["baseline", "hybrid"]
			hyper_parameter_search(data_CEL, labelled_data, k, interesting="hybrid", N=50, metric="FalsePositiveRate(levels=$levels)", max=false)
		catch
			continue
		end
	end

	train_tss = vcat(Bchoices, Hchoices)
	n = length(train_tss)
	ind = [data_CEL.shot_dict[ts] for ts in train_tss]

	dat = data_CEL.cosine_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

	cos_mode = mode(CD)
	cos_ext = extrema(CD)

	dat = data_CEL.flat_top_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

	ft_mode = mode(CD)
	ft_ext = extrema(CD)

	shots = [j for (i, j) in res.ts]
	labels = (@subset all_2D_features @byrow begin
		:tok == "aug"
		in(:shot, shots)
	end).current_heating

	BACC = 0
	for S in 1:50	
		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		labelled_ts = collect(labelled_data |> keys)
		labelled_y = collect(labelled_data |> values)
		shot_dict = Dict([a => n for (n, a) in enumerate(data_CEL.tok_shots)])
		labelled_ind = [shot_dict[k] for k in labelled_ts]

		shot_dict = Dict([a => n for (n, a) in enumerate(data_CEL.tok_shots)])
		labels = collect(labelled_data |> values) |> unique

		train_ind = training_partion(labelled_data, labels, k, S=S)

		X = exp.(-(Array(data_CEL.cosine_cost[!, 2:end]) ./ cos_med).^2) .* 
			exp.(-(Array(data_CEL.flat_top_cost[!, 2:end]) ./ ft_med).^2)

		K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
		model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed)

		KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
		ỹ, _ = svmpredict(model, KK)
		# ỹ, _ = KNN(KK, labelled_y[train_ind])
		BACC += BalancedAccuracy()(ỹ, labelled_y[Not(train_ind)])
	end
	BACC = BACC / 50

	append!(meta_D_BH, 
		DataFrame(:parameters => para,
			:metric => ACC,
			:Bacc => BACC,
			:COS => cos_med,
			:cos_hist_mode => cos_mode,
			:cos_hist_ext => cos_ext,
			:FT => ft_med,
			:ft_hist_mode => ft_mode,
			:ft_hist_ext => ft_ext
			))
end
meta_D_BH |> vscodedisplay
# data_BH = let
#     BH_features = ["IP", "PNBI", "BETAPOL", "Q95"]
#     BH_tss = tok_shots((@subset which(BH_features) @byrow in(:tok, ["aug"])))
#     DTW_hyp_1(BH_tss, BH_features, 10; L=100)
# end

wanted_para = "BETAPOL_IP_LI_NGW_PECRH_PICRH_PNBI_Q95"
COS_med = 59.47
FT_med = 1429.65
data_BH = data_bhind[wanted_para]
# RES, CONFIDENCE, KERNEL, (COS_med, FT_med)= let
# 	labelled_data = let
# 	    labelled_shots = vcat(Bchoices, Hchoices)
# 	    labels = vcat(repeat(["baseline"], length(Bchoices)),
# 	        repeat(["hybrid"], length(Hchoices))) 
# 	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
# 	end
# 	levels = ["baseline", "hybrid"]
# 	hyper_parameter_search(data_bhind["BETAPOL_IP_LI_NGW_PECRH_PICRH_PNBI_Q95"], labelled_data, [10, 10], interesting="hybrid", N=100, metric="FalsePositiveRate(;levels=$levels)", max=false)
# end
RES, CONFIDENCE, KERNEL, MODEL = let
	labelled_data = let
	    labelled_shots = vcat(Bchoices, Hchoices)
	    labels = vcat(repeat(["baseline"], length(Bchoices)),
	        repeat(["hybrid"], length(Hchoices))) 
	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	classify!(data_bhind[wanted_para], labelled_data, (COS_med, FT_med), interesting="hybrid")
end

# append!(RES, DataFrame(:ts => vcat(Bchoices, Hchoices), :predict => vcat(["baseline" for _ in 1:length(Bchoices)], ["hybrid" for _ in 1:length(Hchoices)])))
# RES.ITPA .= ""
# AUG_W = (@subset SELDB5.original_space.ELMy @byrow in(:TOK, ["AUG", "AUGW"]))
# gbdf, keyz = grouping(AUG_W, :SHOT)

# shot_dict = Dict([j for (i, j) in RES.ts] .=> 1:size(RES, 1))
# for (df, key) in zip(gbdf, keyz)
# 	if !in(df.SHOT[1], AUG_shots)
# 		continue
# 	elseif in("YES", df.HYBRID)
# 		RES[shot_dict[key[1]], :ITPA] = "YES"
# 	else
# 		RES[shot_dict[key[1]], :ITPA] = "NO"
# 	end
# end

# RES
# @subset RES @byrow begin
# 	:predict == "hybrid"
# 	:ITPA == "YES"
# end

let
	k = [10, 10]

	train_tss = vcat(Bchoices, Hchoices)
	n = length(train_tss)

	labelled_data = let
	    labelled_shots = vcat(Bchoices, Hchoices)
	    labels = vcat(repeat(["baseline"], length(Bchoices)),
	        repeat(["hybrid"], length(Hchoices))) 
	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

	train_ind = training_partion(labelled_data, unique(labelled_y), k, S=10)

	dat = Array(data_BH.cosine_cost[labelled_ind[train_ind], 2:end])
	mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ COS_med).^2)

	dat = Array(data_BH.flat_top_cost[labelled_ind[train_ind], 2:end])
	mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ FT_med).^2)

	mat = mat_c .* mat_ft

	f = Figure();
	a1 = Axis(f[1, 1])
	# a1 = Axis(f[1:3, 1:20], title="", ylabel="", ylabelsize=20, yticks=1:50:500, ylabelpadding=20)
	# a2 = Axis(f[1:3, 21:25], yticklabelsvisible=false, yticksvisible=false, title="")
	he = heatmap!(a1, mat, colorscale=exp10, colormap=:oslo)

	cntmap = countmap(labelled_y[Not(train_ind)])
	lines!(a1, [k[1]+0.5], [0.5, cntmap["baseline"]+cntmap["hybrid"]+0.5], color=:white, linewidth=4)
	# lines!(a1, [k[1]+k[2]+0.5], [cntmap["C-LH"], length(labelled_y[Not(train_ind)])].+0.5, color=:white, linewidth=4)

	lines!([0, k[1]+k[2]].+0.5, [cntmap["baseline"]+0.5], color=:white, linewidth=4)
	# lines!([k[1], length(labelled_y[train_ind])].+0.5, [cntmap["C-LH"]+cntmap["C-EH"]+0.5], color=:white, linewidth=4)

	# lines!([k[1], k[1]].+0.5, [0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)
	# lines!([k[1]+k[2], k[1]+k[2]].+0.5, [k[1]+0.5, +(k...)+0.5], color=:white, linewidth=4)
	# lines!([0.5, k[1]+k[2]+0.5], [k[1], k[1]].+0.5, color=:white, linewidth=4)
	# lines!([k[1]+0.5, +(k...)+0.5], [k[1]+k[2]+0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)

	Colorbar(f[1, 2], he)
	# a1.xticks = (1:length(Bchoices), ["#$j" for (i,j) in Bchoices])
	a1.yticks = (1:length(labelled_ts[Not(train_ind)]), ["#$j" for (i, j) in labelled_ts[Not(train_ind)]])
	a1.xticks = (1:length(labelled_ts[train_ind]), ["#$j" for (i, j) in labelled_ts[train_ind]])
	a1.xticklabelrotation = π/5
	
	# a1.xticklabelrotation = π/2
	# a2.xticklabelrotation = π/2
	
	# a1.xticklabelsize = 15
	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/BH_kernel_matrix.png", f)
	f
end
let
	# CairoMakie.activate!()
	k = [10, 10]
	
	train_tss = vcat(Bchoices, Hchoices)
	n = length(train_tss)

	labelled_data = let
	    labelled_shots = vcat(Bchoices, Hchoices)
	    labels = vcat(repeat(["baseline"], length(Bchoices)),
	        repeat(["hybrid"], length(Hchoices))) 
	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	labelled_ts = collect(labelled_data |> keys)
    labelled_y = collect(labelled_data |> values)
    shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])
    labelled_ind = [shot_dict[k] for k in labelled_ts]

	train_ind = training_partion(labelled_data, unique(labelled_y), k)
	n = sum(train_ind)

	fig = Figure();
	ax = Axis(fig[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20)
	ax1 = Axis(fig[1, 2], xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20)
	ax2 = Axis(fig[1, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20)

	dat = Array(data_BH.cosine_cost[labelled_ind[train_ind], 2:end])
	mat_c = exp.(-(dat[:, labelled_ind[train_ind]] ./ COS_med).^2)
	CD = vcat([mat_c[i, i+1:end] for i in 1:n-1]...)
	hist!(ax, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
	# hlines!(ax, [COS_med], xmin=0.05, xmax=0.95, color=:red)

	dat = Array(data_BH.flat_top_cost[labelled_ind[train_ind], 2:end])
	mat_ft = exp.(-(dat[:, labelled_ind[train_ind]] ./ FT_med).^2)
	CD = vcat([mat_ft[i, i+1:n] for i in 1:n-1]...)
	hist!(ax1, CD, normalization=:probability, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
	# # hlines!(ax1, [FT_med], xmin=0.05, xmax=0.95, color=:red)
	
	mat = mat_c .* mat_ft
	CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
	hist!(ax2, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
	# # # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/cost_spread_BH.png", fig)
	fig
end
CairoMakie.activate!()
let
	pred_shots = [j for (i, j) in RES.ts]
	shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])

	CH = [current_heating("aug", shot) for shot in pred_shots]
	
	col_dict = Dict("baseline"=>1, "hybrid"=>2, "training"=>3)
	CO_dict = Dict("CO"=> :cross, "C-EH"=> :ltriangle, "C-LH"=> :rtriangle)

	training_D = DataFrame(:ts => vcat(Bchoices, Hchoices),
						:predict => vcat(["training" for _ in 1:length(Bchoices)],
							["training" for _ in 1:length(Hchoices)]),
						:conf => vcat([(ts[2], (0.6*cos(n+1))) for (n, ts) in enumerate(Bchoices)],
									[(ts[2], (0.4*cos(n+1))) for (n, ts) in enumerate(Hchoices)]))
	training_D.CH = [current_heating(ts...) for ts in training_D.ts]

	RES.conf = [Tuple.(zip(ts[2], CONFIDENCE[1, n])) for (n, ts) in enumerate(RES.ts)]
	
	# RES_updated = vcat(RES[!, [:ts, :predict, :conf]], training_D)
	RES.CH .= CH
	ỹ = RES.predict
	
	HYBRID = (@subset RES @byrow :predict == "hybrid")[1:4:end, :]
	CO = (@subset RES @byrow :CH == "CO")
	
	f = Figure(size=(1200, 700));
	a0 = Axis(f[1:3, 1:12])
	a = Axis(f[4:16, 1:12])
	a1 = Axis(f[17:19, 1:12])
    a_sp = scatter!(a, RES.conf; colormap=:tab10, 
		color=[col_dict[el] for el in RES.predict], 
		colorrange = (1, 3),
		strokewidth = 0.1, 
		markersize = 18,
		marker = [CO_dict[label] for label in RES.CH],
		label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ])
	a0.xticklabelsvisible = false
	a0.yticklabelsvisible = false
	hidespines!(a0)
	a0.limits = (nothing, (-0.7, 0.7))
	hidedecorations!(a0)
	a.xticks=(10000:10000:40000, ["10,000", "20,000", "30,000", "40,000"])
	a1.yticklabelsvisible = false
	hidespines!(a1)
	a1.limits = (nothing, (-0.6, 0.6))
	hidedecorations!(a1)
	
	linkxaxes!(a0, a, a1)
	
	a0_sp = scatter!(a0, training_D.conf[1:length(Bchoices)]; colormap=:tab10, 
		color=3, 
		colorrange = (1, 3),
		strokewidth = 0.1, 
		markersize = 18,
		marker = [CO_dict[label] for label in training_D.CH[1:length(Bchoices)]],
		label = "training")

	scatter!(a1, training_D.conf[length(Bchoices)+1:end]; colormap=:tab10, 
		color=3, 
		colorrange = (1, 3),
		strokewidth = 0.1, 
		markersize = 18,
		marker = [CO_dict[label] for label in training_D.CH[length(Bchoices)+1:end]],
		label = "training")

    pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
    for (n, (ts, coord)) in enumerate(zip(HYBRID.ts, HYBRID.conf))
        text!(a, coord, text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13)
    end
	pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
    for (n, (ts, coord)) in enumerate(zip(CO.ts, CO.conf))
        text!(a, coord, text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=13)
    end
	pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
	tr_D_CO = @subset training_D @byrow :CH == "CO"
    for (n, (ts, coord)) in enumerate(zip(tr_D_CO.ts, tr_D_CO.conf))
		text!(a1, coord, text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13)
    end

	# extra = (@subset training_data @byrow in(:ts, [("aug", shot) for shot in [12000, 12032]]))
	# for (n, (ts, coord)) in enumerate(zip(extra.ts, extra.conf))
	#        text!(a, coord, text="$(ts[2])", align=pos[mod(n,2+2)+1], fontsize=7)
	#    end

	# a.xgridvisible = false
	# a.ygridvisible = false
	
	group_marker = [MarkerElement(marker = marker, color = :black,
		strokecolor = :transparent,
		markersize = 12) for marker in [:ltriangle, :rtriangle, :cross]]

	Legend(f[15, 2],
		group_marker,
		["EH", "LH", "CO"],
		framevisible=false, labelsize=20)

	colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
	group_marker = [MarkerElement(marker = :circ, color = colors_3[i],
	    strokecolor = :transparent,
	    markersize = 12) for i in 1:3]

    Legend(f[15, 1],
	    group_marker,
	    ["baseline", "hybrid", "training"],
	 	framevisible=false, labelsize=20)
	# axislegend(a0, unique=true, merge=true, position=:lb, labelsize=20)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/Baseline_Hybrid_classified.png", f)
	
	println("hybrid = ", (@subset RES @byrow :predict == "hybrid").CH |> countmap)
	println("baseline = ", (@subset RES @byrow :predict == "baseline").CH |> countmap)
	
    f
end

figs = Figure[]
for (wanted_para, metric, COS_med, FT_med) in eachrow(meta_D_BH[!, [:parameters, :metric, :COS, :FT]])
	RES, CONFIDENCE, KERNEL, MODEL = let
		labelled_data = let
			labelled_shots = vcat(Bchoices, Hchoices)
			labels = vcat(repeat(["baseline"], length(Bchoices)),
				repeat(["hybrid"], length(Hchoices))) 
			d = OrderedDict([ts for ts in labelled_shots] .=> labels)
		end
		classify!(data_bhind[wanted_para], labelled_data, (COS_med, FT_med), interesting="hybrid")
	end
	f = let 
		f = Figure(size=(1000, 1200));
		gl1 = f[1, 1] = GridLayout()
		# gl1[1, 1] = Axis(f)
		try let
				k = [10, 10]
		
				train_tss = vcat(Bchoices, Hchoices)
				n = length(train_tss)
		
				labelled_data = let
					labelled_shots = vcat(Bchoices, Hchoices)
					labels = vcat(repeat(["baseline"], length(Bchoices)),
						repeat(["hybrid"], length(Hchoices))) 
					d = OrderedDict([ts for ts in labelled_shots] .=> labels)
				end
				labelled_ts = collect(labelled_data |> keys)
				labelled_y = collect(labelled_data |> values)
				shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])
				labelled_ind = [shot_dict[k] for k in labelled_ts]
		
				train_ind = training_partion(labelled_data, unique(labelled_y), k, S=10)
		
				dat = Array(data_BH.cosine_cost[labelled_ind[train_ind], 2:end])
				mat_c = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ COS_med).^2)
		
				dat = Array(data_BH.flat_top_cost[labelled_ind[train_ind], 2:end])
				mat_ft = exp.(-(dat[:, labelled_ind[Not(train_ind)]] ./ FT_med).^2)
		
				mat = mat_c .* mat_ft
				
				# a1 = Axis(f[1:3, 1:20], title="", ylabel="", ylabelsize=20, yticks=1:50:500, ylabelpadding=20)
				# a2 = Axis(f[1:3, 21:25], yticklabelsvisible=false, yticksvisible=false, title="")
				ax = Axis(gl1[1, 1])
				he = heatmap!(ax, mat, colorscale=exp10, colormap=:oslo)
		
				cntmap = countmap(labelled_y[Not(train_ind)])
				lines!(ax, [k[1]+0.5], [0.5, cntmap["baseline"]+cntmap["hybrid"]+0.5], color=:white, linewidth=4)
				# lines!(a1, [k[1]+k[2]+0.5], [cntmap["C-LH"], length(labelled_y[Not(train_ind)])].+0.5, color=:white, linewidth=4)
		
				lines!(ax, [0, k[1]+k[2]].+0.5, [cntmap["baseline"]+0.5], color=:white, linewidth=4)
				# lines!([k[1], length(labelled_y[train_ind])].+0.5, [cntmap["C-LH"]+cntmap["C-EH"]+0.5], color=:white, linewidth=4)
		
				# lines!([k[1], k[1]].+0.5, [0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)
				# lines!([k[1]+k[2], k[1]+k[2]].+0.5, [k[1]+0.5, +(k...)+0.5], color=:white, linewidth=4)
				# lines!([0.5, k[1]+k[2]+0.5], [k[1], k[1]].+0.5, color=:white, linewidth=4)
				# lines!([k[1]+0.5, +(k...)+0.5], [k[1]+k[2]+0.5, k[1]+k[2]+0.5], color=:white, linewidth=4)
		
				Colorbar(gl1[1, 2], he)
				# a1.xticks = (1:length(Bchoices), ["#$j" for (i,j) in Bchoices])
				ax.yticks = (1:length(labelled_ts[Not(train_ind)]), ["#$j" for (i, j) in labelled_ts[Not(train_ind)]])
				ax.xticks = (1:length(labelled_ts[train_ind]), ["#$j" for (i, j) in labelled_ts[train_ind]])
				ax.xticklabelrotation = π/5
				
				# a1.xticklabelrotation = π/2
				# a2.xticklabelrotation = π/2
				
				# a1.xticklabelsize = 15
				# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/BH_kernel_matrix.png", f)
			end
		catch
			continue
		end
		
		gl2 = GridLayout()
		gl2[1, 1:3] = [Axis(f, xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20),
			Axis(f, xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20),
			Axis(f, xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20)]
		try let
				# CairoMakie.activate!()
				k = [10, 10]
				
				train_tss = vcat(Bchoices, Hchoices)
				n = length(train_tss)
		
				labelled_data = let
					labelled_shots = vcat(Bchoices, Hchoices)
					labels = vcat(repeat(["baseline"], length(Bchoices)),
						repeat(["hybrid"], length(Hchoices))) 
					d = OrderedDict([ts for ts in labelled_shots] .=> labels)
				end
				labelled_ts = collect(labelled_data |> keys)
				labelled_y = collect(labelled_data |> values)
				shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])
				labelled_ind = [shot_dict[k] for k in labelled_ts]
		
				train_ind = training_partion(labelled_data, unique(labelled_y), k)
				n = sum(train_ind)
		
				# fig = Figure();
				# ax2 = Axis(f[1, 2])
				# ax2 = Axis(fig[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20)
				# ax1 = Axis(fig[1, 2], xticks=(1.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20)
				# ax2 = Axis(fig[1, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20)
		
				dat = Array(data_BH.cosine_cost[labelled_ind[train_ind], 2:end])
				mat_c = exp.(-(dat[:, labelled_ind[train_ind]] ./ COS_med).^2)
				CD = vcat([mat_c[i, i+1:end] for i in 1:n-1]...)
				hist!(gl2[1, 1], CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
				# hlines!(ax, [COS_med], xmin=0.05, xmax=0.95, color=:red)
		
				dat = Array(data_BH.flat_top_cost[labelled_ind[train_ind], 2:end])
				mat_ft = exp.(-(dat[:, labelled_ind[train_ind]] ./ FT_med).^2)
				CD = vcat([mat_ft[i, i+1:n] for i in 1:n-1]...)
				hist!(gl2[1, 2], CD, normalization=:probability, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
				# # hlines!(ax1, [FT_med], xmin=0.05, xmax=0.95, color=:red)
				
				mat = mat_c .* mat_ft
				CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
				hist!(gl2[1, 3], CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray30)
				# # # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/cost_spread_BH.png", fig)
				# fig
			end
		catch
			continue
		end
		
		gl3 = f[2:3, 1:2] = GridLayout()
		# gl3[1:3, 1] = [Axis(f) for _ in 1:3]
		# rowsize!(gl3, 1, Relative(0.1))
		# rowsize!(gl3, 3, Relative(0.1))
		try let
				pred_shots = [j for (i, j) in RES.ts]
				shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])
		
				CH = [current_heating("aug", shot) for shot in pred_shots]
				
				col_dict = Dict("baseline"=>1, "hybrid"=>2, "training"=>3)
				CO_dict = Dict("CO"=> :cross, "C-EH"=> :ltriangle, "C-LH"=> :rtriangle)
		
				training_D = DataFrame(:ts => vcat(Bchoices, Hchoices),
									:predict => vcat(["training_baseline" for _ in 1:length(Bchoices)],
										["training_hybrid" for _ in 1:length(Hchoices)]),
									:conf => vcat([(ts[2], (0.6*cos(n+1))) for (n, ts) in enumerate(Bchoices)],
												[(ts[2], (0.55*cos(n+1))) for (n, ts) in enumerate(Hchoices)]))
				training_D.CH = [current_heating(ts...) for ts in training_D.ts]
		
				RES.conf = [Tuple.(zip(ts[2], CONFIDENCE[1, n])) for (n, ts) in enumerate(RES.ts)]
				
				# RES_updated = vcat(RES[!, [:ts, :predict, :conf]], training_D)
				RES.CH .= CH
				ỹ = RES.predict
				
				HYBRID = (@subset RES @byrow :predict == "hybrid")[1:2:end, :]
				CO = (@subset RES @byrow begin
					:CH == "CO"
					:predict == "baseline"
				end)
		
				axtop = Axis(gl3[1:3, 1:12])
				axmid = Axis(gl3[4:16, 1:12], ylabel="confidence", xlabel="shot #")
				axbot = Axis(gl3[17:19, 1:12])
				a_sp = scatter!(axmid, RES.conf; colormap=:tab10, 
					color=[col_dict[el] for el in RES.predict], 
					colorrange = (1, 3),
					strokewidth = 0.1, 
					markersize = 18,
					marker = [CO_dict[label] for label in RES.CH],
					label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ])
				axtop.xticklabelsvisible = false
				axtop.yticklabelsvisible = false
				hidespines!(axtop)
				axtop.limits = (nothing, (-0.9, 0.9))
				hidedecorations!(axtop)
				axmid.xticks=(10000:10000:40000, ["10,000", "20,000", "30,000", "40,000"])
				axbot.yticklabelsvisible = false
				hidespines!(axbot)
				axbot.limits = (nothing, (-0.9, 0.9))
				hidedecorations!(axbot)
				linkxaxes!(axtop, axmid, axbot)
				
				a0_sp = scatter!(axtop, training_D.conf[1:length(Bchoices)]; colormap=:tab10, 
					color=3, 
					colorrange = (1, 3),
					strokewidth = 0.1, 
					markersize = 18,
					marker = [CO_dict[label] for label in training_D.CH[1:length(Bchoices)]],
					label = "training")
		
				scatter!(axbot, training_D.conf[length(Bchoices)+1:end]; colormap=:tab10, 
					color=3, 
					colorrange = (1, 3),
					strokewidth = 0.1, 
					markersize = 18,
					marker = [CO_dict[label] for label in training_D.CH[length(Bchoices)+1:end]],
					label = "training")
		
				pos = [(:left, :top), (:right, :bottom), (:left, :bottom), (:right, :top)]
				for (n, (ts, coord)) in enumerate(zip(HYBRID.ts, HYBRID.conf))
					text!(axmid, coord, text="$(ts[2])", align=pos[mod(n+0,4)+1], fontsize=13)
				end
				pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
				for (n, (ts, coord)) in enumerate(zip(CO.ts, CO.conf))
					text!(axmid, coord, text="$(ts[2])", align=pos[mod(n+0,4)+1], fontsize=13)
				end
				pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
				tr_D_hy = @subset training_D @byrow :predict == "training_hybrid"
				for (n, (ts, coord)) in enumerate(zip(tr_D_hy.ts, tr_D_hy.conf))
					text!(axbot, coord, text="$(ts[2])", align=pos[mod(n+0,4)+1], fontsize=13)
				end
				pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
				tr_D_CO = @subset training_D @byrow :predict == "training_baseline"
				for (n, (ts, coord)) in enumerate(zip(tr_D_CO.ts, tr_D_CO.conf))
					text!(axtop, coord, text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13)
				end
		
				group_marker = [MarkerElement(marker = marker, color = :black,
				strokecolor = :transparent,
				markersize = 12) for marker in [:ltriangle, :rtriangle, :cross]]
		
				Legend(f[2:3, 1:2][14, 1],
					group_marker,
					["EH", "LH", "CO"],
					framevisible=false, labelsize=20, halign=:left, valign=:top)
		
				colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
				group_marker = [MarkerElement(marker = :circ, color = colors_3[i],
					strokecolor = :transparent,
					markersize = 12) for i in 1:3]
		
				Legend(f[2:3, 1:2][17, 1],
					group_marker,
					["baseline", "hybrid", "training"],
					framevisible=false, labelsize=20, halign=:left, valign=:top)
		
				# # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/Baseline_Hybrid_classified.png", f)
				
				# # println("hybrid = ", (@subset RES @byrow :predict == "hybrid").CH |> countmap)
				# # println("baseline = ", (@subset RES @byrow :predict == "baseline").CH |> countmap)
				
				# f
			end
		catch
			continue
		end
		
		# f.layout[1, 1] = gl1
		f.layout[1, 2] = gl2
		# f.layout[2:3, 1:2] = gl3
		Label(f[0, :], text = *(split(wanted_para, "_"), delim=" "), fontsize = 20)
		
		f
	end
	push!(figs, f)
end
# save("/Users/joe/Project/PhD/EuroFusion/EUROfusion_ML_2024/Presentations/Multidimensional_qunatity_classification_20_01_25/Baseline_Hybrid_classified.png", figs[end])
figs[6]


# for (CLn, CEn, COn) in zip([60, 35, 10], [30, 20, 10], [5, 10, 10])
# 	for S in [14, 123, 69]
# 		CLchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false))
# 		CEchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false))
# 		COchoices = vcat(sample(Random.seed!(S), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false))

# 		k = floor.(Int, 0.6.*[CLn, CEn, COn])
# 		res, confidence, kernel_matrix, (cos_med, ft_med), ACC = let
# 			labelled_data = let
# 				labelled_shots = vcat(CLchoices, CEchoices, COchoices)
# 				labels = vcat(repeat(["C-LH"], length(CLchoices)),
# 					repeat(["C-EH"], length(CEchoices)),
# 					repeat(["CO"], length(COchoices))) 
# 				d = OrderedDict([ts for ts in labelled_shots] .=> labels)
# 			end
# 			hyper_parameter_search(data_CEL, labelled_data, k, interesting="CO", N=60)
# 		end

# 		train_tss = vcat(CLchoices, CEchoices, COchoices)
# 		n = length(train_tss)
# 		ind = [data_CEL.shot_dict[ts] for ts in train_tss]

# 		dat = data_CEL.cosine_cost
# 		mat = Array(dat[ind, ind.+1])
# 		CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

# 		cos_mode = mode(CD)
# 		cos_ext = extrema(CD)

# 		dat = data_CEL.flat_top_cost
# 		mat = Array(dat[ind, ind.+1])
# 		CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))	

# 		ft_mode = mode(CD)
# 		ft_ext = extrema(CD)

# 		shots = [j for (i, j) in res.ts]
# 		labels = (@subset all_2D_features @byrow begin
# 			:tok == "aug"
# 			in(:shot, shots)
# 		end).current_heating

# 		acc, b_acc = Accuracy()(res.predict, labels), BalancedAccuracy()(res.predict, labels)

# 		append!(meta_D, 
# 			DataFrame(:parameters => "IP_PNBI", 
# 				:random_seed => S,
# 				:LH => CLn,
# 				:EH => CEn,
# 				:CO => COn,
# 				:Bacc_cv => ACC,
# 				:acc_D => acc,
# 				:Bacc_D => b_acc,
# 				:COS => cos_med,
# 				:cos_hist_mode => cos_mode,
# 				:cos_hist_ext => cos_ext,
# 				:FT => ft_med,
# 				:ft_hist_mode => ft_mode,
# 				:ft_hist_ext => ft_ext
# 				))
# 	end
# end
COS_med, FT_med = 20.3 ,337.7
RES, CONFIDENCE, KERNEL = let
	labelled_dict = let
		labelled_shots = vcat(Bchoices, Hchoices)
		labels = vcat(repeat(["baseline"], length(Bchoices)),
			repeat(["hybrid"], length(Hchoices))) 
		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	classify!(data_BH, labelled_dict, (COS_med, FT_med), interesting="hybrid")
end

shots = [j for (i, j) in data_BH.tok_shots]

res = vcat(DataFrame(:ts => vcat(Bchoices, Hchoices), :predict => vcat(repeat(["baseline"], Bn), repeat(["hybrid"], Hn))), RES)
label_dict = Dict([j for (i, j) in res.ts] .=> res.predict)


D = deepcopy(@subset original_space(DATA_) @byrow begin
	in(:TOK, ["AUG", "AUGW"])
	in(:SHOT, shots)
end)
D.BH .= ""
gbdf, keyz = grouping(D, :SHOT)
for key in keyz
	gbdf[key].BH .= label_dict[key[1]]
end
gbdf
D = vcat(gbdf...)
using StatsPlots
corrplot(Array(D[!, [:IP, :BT, :NEL, :PLTH]]))

dict = Dict("NO"=>"baseline", "YES"=>"hybrid", "HYBRID"=>"hybrid")
D.HYBRID = [dict[el] for el in D.HYBRID]

gbdf, keyz = grouping(D, [:BH, :HYBRID])
keyz[4]
gbdf[keyz[4]]



(@subset D @byrow :BH == "baseline")
reg = Regression((@subset D @byrow :BH == "baseline"), ols(), single_machine())
reg.results
reg.results

# features = ["IP", "BETAPOL", "Q95", "NGW", "PNBI", "PECRH", "PICRH"]

# tss = tok_shots((@subset which(features) @byrow in(:HYBRID_v1, ["YES", "NO"])))
# data = DTW_hyp_1(tss, features, 10, L=100)

# 1NN function
function KNN(kernel::Array, train_labels::Vector{String})
    tr_ℓ, te_ℓ = size(kernel)
    @assert length(train_labels) == tr_ℓ

    dict = Dict(1:tr_ℓ .=> train_labels)

    max_val = Vector{Number}(undef, te_ℓ)
    labels = Vector{Int}(undef, te_ℓ)
    for n in 1:te_ℓ
        max_val[n], labels[n] = findmax(kernel[:, n])
    end
    return [dict[element] for element in labels], max_val
end

features = ["IP", "PNBI"]
tss = [("aug", shot) for shot in AUG_shots]
data = DTW_hyp_1(tss, features, 10, L=100)

Random.seed!(17)
CLn = 20
CEn = 20
COn = 5
CLchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false)
CEchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false)
COchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false)

CEL_feat = (@subset all_2D_features @byrow begin
    in(:tok, [i for (i, j) in tss])
    in(:shot, [j for (i, j) in tss])
end)
cardinality_metadata(CEL_feat, :current_heating) 


labelled_data = let
    labelled_shots = vcat(CLchoices, CEchoices, COchoices)
    labels = vcat(repeat(["C-LH"], length(CLchoices)),
        repeat(["C-EH"], length(CEchoices)),
        repeat(["CO"], length(COchoices))) 
    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
end

res, confidence, kernel_matrix = hyper_parameter_search(data, labelled_data, [9,9,3], interesting="CO", N=5)

shots = [j for (i, j) in res.ts]
labels = (@subset all_2D_features @byrow begin
    :tok == "aug"
    in(:shot, shots)
end).current_heating
res.label .= labels
res.conf = Tuple.(zip(confidence[1, :], confidence[2, :]))
curr_over = (@subset res @byrow :predict == "CO")

acc, b_acc = Accuracy()(labels, res.predict), BalancedAccuracy()(labels, res.predict)
println("accuracy = $(acc), balanced accuracy = $(b_acc)")

col_dict = Dict("CO"=>1, "C-EH"=>2, "C-LH"=>3)
ỹ = res.predict
y = res.label

f, a, s = scatter(confidence[1, :], confidence[2, :]; colormap=:tab10, 
    color=[col_dict[el] for el in ỹ],
    strokewidth = 0.1,
    label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
    markersize = [(i == j ? 20 : 10) for (i, j) in zip(y, ỹ)],
    marker = [(i == j ? '∘' : (:utriangle)) for (i, j) in zip(y, ỹ)]
)

pos = [(:right, :bottom), (:left, :bottom), (:right, :top), (:left, :top)]
for (n, (ts, coord)) in enumerate(zip(curr_over.ts, curr_over.conf))
    text!(a, coord, text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=7)
end

axislegend(a, unique=true, merge=true, position=:lt)

# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/CO_EH_LH_classification.png", f)
f









using CairoMakie


let
    CairoMakie.activate!()
	shot = 13744
	x1_lim = 8.6
	x2_lim = 6.5
	fig1 = Figure(size=(1900, 900))
	ax1 = Axis(fig1[1:3, 2:4], title="#$shot", titlesize=25)
	ax2 = Axis(fig1[4:6, 2:4])

	features = ["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"]
	t = FTIP[("aug", shot)]
    
	for (n, (feat, norm, axis, label)) in enumerate(zip(features, [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		BV = t[1] .< P.t .< t[2]
		Y = abs.(P.y.y .* norm)
		x_shade = P.t[BV]
		y_shade = Y[BV]
		band!(axis, x_shade, 0, y_shade, color=(:cyan, 0.15))
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=:dracula, colorrange=(1, 8), color=n, label=label)
	end
	ax1.limits = ((0, x1_lim), (-0.1, 2.5))
	ax2.limits = ((0, x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)
	
	shot = 8045
	t = FTIP[("aug", shot)]
	ax3 = Axis(fig1[1:3, 5:7], title="#$shot", titlesize=25)
	ax4 = Axis(fig1[4:6, 5:7])
	for (n, (feat, norm, axis, label)) in enumerate(zip(features, [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		BV = t[1] .< P.t .< t[2]
		Y = abs.(P.y.y .* norm)
		x_shade = P.t[BV]
		y_shade = Y[BV]
		band!(axis, x_shade, 0, y_shade, color=(:cyan, 0.15))
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=:dracula, colorrange=(1, 8), color=n, label=label)
	end
	ax3.limits = ((0, x2_lim), (-0.1, 2.5))
	ax4.limits = ((0, x2_lim), (-0.1, 2.5))
	ax3.yticks = 0:0.5:2.5
	ax4.yticks = 0:0.5:2.5
	# fig1[1:2, 8] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	# fig1[4:5, 8] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/shot_comparison_stochastic_FTNBI.png", fig1)
	fig1
end

@subset original_space(DATA_) @byrow begin
    :TOK == "AUG"
    :SHOT == 27910
end

all_IB = [8127, 29619, 29620, 29621, 29622, 29623, 29636, 29953, 29956, 29957, 29958, 29962, 29963, 
    29964, 29965, 29976, 29977, 29979, 29980, 29981, 32103, 32104, 32105, 32106, 32109, 32110, 32111, 32112, 32113, 32114, 32115, 
    32120, 32121, 32122, 32123, 32463, 32464, 32471, 32472, 32474, 32475, 32489, 33370, 33371, 33374, 33375, 33376, 33377, 33406, 
    33407, 33409, 33427, 34441, 34442, 34443, 34444, 34445, 34446, 34449, 34450, 34454, 34489, 34490, 34838, 34839, 34840, 34841, 
    34842, 34844, 35975, 36176, 36177, 36178, 37080, 37081, 37082, 37093, 37094, 37096, 39124, 39125, 39126, 39128, 39129, 40202, 
    40203, 40206, 40207, 40410, 40411, 40851
]

Bchoices = [("aug", shot) for shot in all_IB]
Hchoices = [("aug", shot) for shot in [11190, 16688, 16736, 18046, 18869, 18880, 19314, 25764, 26338, 26913, 27930, 34769, 34770, 
    34774, 34775, 34776, 36443, 36408]]


Bn = length(Bchoices)
Hn = length(Hchoices)

d = DataFrame(:ts => vcat(Bchoices, Hchoices), :label => vcat(repeat(["baseline"], Bn), repeat(["hybrid"], Hn)))
cardinality_metadata(d, :label) |> clipboard


RES, CONFIDENCE, KERNEL, (CD, MD) = let
	labelled_dict = let
	    labelled_shots = vcat(Bchoices, Hchoices)
	    labels = vcat(repeat(["baseline"], length(Bchoices)),
	        repeat(["hybrid"], length(Hchoices))) 
	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	hyper_parameter_search(data_main, labelled_dict, [35, 10], interesting="hybrid", N=200)
end
let
	shots = [j for (i, j) in RES.ts]
	CH = (@subset all_2D_features @byrow begin
	    :tok == "aug"
	    in(:shot, shots)
	end).current_heating
	
	col_dict = Dict("baseline"=>1, "hybrid"=>2)
	CO_dict = Dict("CO"=>1, "C-EH"=>2, "C-LH"=>3)
	
    ỹ = RES.predict
	RES.conf = Tuple.(zip(1:size(CONFIDENCE, 2), CONFIDENCE[1, :]))
	RES.CH .= CH
	
	HYBRID = (@subset RES @byrow :predict == "hybrid")
    f, a, s = scatter(1:size(CONFIDENCE, 2), CONFIDENCE[1, :]; colormap=:tab10, 
        color=[col_dict[el] for el in ỹ],
		strokewidth = 0.1,
        label = [label => (;colormap=:tab10, colorrange=(1, 2), color=col_dict[label]) for label in ỹ])
	
    pos = [(:right, :bottom), (:left, :bottom), (:right, :top), (:left, :top)]
    for (n, (ts, coord)) in enumerate(zip(HYBRID.ts, HYBRID.conf))
        text!(a, coord, text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=7)
    end

	# extra = (@subset RES @byrow in(:ts, [("aug", shot) for shot in [12000, 12032]]))
	# for (n, (ts, coord)) in enumerate(zip(extra.ts, extra.conf))
 #        text!(a, coord, text="$(ts[2])", align=pos[mod(n,2+2)+1], fontsize=7)
 #    end

    axislegend(a, unique=true, merge=true, position=:lt)
    f
end
let
	train_tss = vcat(Bchoices, Hchoices)
	n = length(train_tss)
	ind = [data_main.shot_dict[ts] for ts in train_tss]

	fig = Figure();
	ax = Axis(fig[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20)
	ax1 = Axis(fig[1, 2], xticks=(0.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{I_p^{80%}}", L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20)

	dat = data_main.cosine_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))
	hist!(ax, CD, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
	# hlines!(ax, [8.85], xmin=0.05, xmax=0.95, color=:red)
    hlines!(ax, [6.32], xmin=0.05, xmax=0.95, color=:red)

	# dat = data_main.magnitude_cost
	# mat = Array(dat[ind, ind.+1])
	# CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))
	# hist!(ax1, CD, scale_to=-0.6, offset=1, direction=:x, color=:gray60)
	# # hlines!(ax1, [1753.64], xmin=0.05, xmax=0.38, color=:red)
	
	dat = data_main.flat_top_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))
	hist!(ax1, CD, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
	# hlines!(ax1, [269.39], xmin=0.6, xmax=0.95, color=:red)
    hlines!(ax1, [134.8], xmin=0.05, xmax=0.95, color=:red)
	
	# # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/cost_spread.png", fig)
	fig
end










Random.seed!(17)
CLn = 20
CEn = 20
COn = 5
CLchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false)
CEchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false)
COchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false)

labelled_data = let
    labelled_shots = vcat(CLchoices, CEchoices, COchoices)
    labels = vcat(repeat(["C-LH"], length(CLchoices)),
        repeat(["C-EH"], length(CEchoices)),
        repeat(["CO"], length(COchoices))) 
    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
end
res, confidence, kernel_matrix = hyper_parameter_search(data_CEL, labelled_data, [10,10,3], interesting="CO", N=50)







function current_heating(tok::String, shot::Int)
	int = (@subset all_2D_features @byrow begin
		:tok == tok
		:shot == shot
	end)
	return int.current_heating[1]
end

all_IB = [7643, 8127, 29619, 29620, 29621, 29622, 29623, 29636, 29953, 29956, 29957, 29958, 29962, 29963, 29964, 29965, 29976, 29977, 29979, 29980, 29981, 32103, 32104, 32105, 32106, 32109, 32110, 32111, 32112, 32113, 32114, 32115, 32120, 32121, 32122, 32123, 32463, 32464, 32471, 32472, 32474, 32475, 32489, 33370, 33371, 33374, 33375, 33376, 33377, 33406, 33407, 33409, 33427, 34441, 34442, 34443, 34444, 34445, 34446, 34449, 34450, 34454, 34489, 34490, 34840, 35975, 36177, 37081, 37094, 37096, 39124, 39125, 39126, 39128, 39129, 40411, 40851, 40410, 37093, 37082, 36178, 37080, 36176, 34842, 40206, 40207, 40202, 40203, 34838, 34839, 34841, 34844
]

chosen_baseline_training = [6483, 7558, 7643, 8007, 8093, 10729, 11126, 13603, 15714, 17389, 20280, 23259, 23971, 
	25719, 25757, 29623, 32225, 33303, 34489, 34490, 35975, 40203, 40411
]

Bchoices = [("aug", shot) for shot in chosen_baseline_training]
Hchoices = [("aug", shot) for shot in [11190, 16688, 16736, 18046, 18869, 18880, 19314, 25764, 26338, 26913, 27930, 34769, 34770, 
	34774, 34775, 34776, 36443, 36408]]

Bn = length(Bchoices)
Hn = length(Hchoices)

data_try = let
    CEL_features = ["IP", "PNBI"]
    CEL_tss = tok_shots((@subset which(CEL_features) @byrow in(:tok, ["aug"])))
    DTW_hyp_1(CEL_tss, CEL_features, 10, L=100)
end

RES, CONFIDENCE, KERNEL, (COS_med, FT_med) = let
	labelled_dict = let
		labelled_shots = vcat(Bchoices, Hchoices)
		labels = vcat(repeat(["baseline"], length(Bchoices)),
			repeat(["hybrid"], length(Hchoices))) 
		d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	hyper_parameter_search(data_try, labelled_dict, [13, 10], interesting="hybrid", N=60)
end

base_lt_0p2 = round(sum(0 .< CONFIDENCE[1, :] .<= 0.2) / sum(0 .< CONFIDENCE[1, :]), digits=2)
base_gt_0p8 = round(sum(0.8 .<= CONFIDENCE[1, :]) / sum(0 .< CONFIDENCE[1, :]), digits=2)
hyb_lt_0p2 = round(sum(-0.2 .<= CONFIDENCE[1, :] .< 0.) / sum(CONFIDENCE[1, :] .< 0), digits=2)
hyb_gt_0p8 = round(sum(CONFIDENCE[1, :] .<= -0.8) / sum(CONFIDENCE[1, :] .< 0), digits=2)
println("COS_med = ", COS_med, "\n", "FT_med = ", FT_med)
println("base %: ($base_lt_0p2, $base_gt_0p8) \n hyb %: ($hyb_lt_0p2, $hyb_gt_0p8)")

let
	train_tss = vcat(Bchoices, Hchoices)
	n = length(train_tss)
	ind = [data_try.shot_dict[ts] for ts in train_tss]

	fig = Figure();
	ax = Axis(fig[1, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20)
	ax1 = Axis(fig[1, 2], xticks=(0.7:1:1.7, [L"\mathrm{DTW_{mag}}|_{I_p^{80%}}", L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20)

	dat = data_try.cosine_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))
	hist!(ax, CD, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
	hlines!(ax, [COS_med], xmin=0.05, xmax=0.95, color=:red)

	dat = data_try.flat_top_cost
	mat = Array(dat[ind, ind.+1])
	CD = abs.(vcat([mat[i, i+1:n] for i in 1:n-1]...))
	hist!(ax1, CD, scale_to=-0.6, offset=2, direction=:x, color=:gray50)
	hlines!(ax1, [FT_med], xmin=0.05, xmax=0.95, color=:red)
	
	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/cost_spread_BH.png", fig)
	fig
end
let
	pred_shots = [j for (i, j) in RES.ts]
	shot_dict = Dict([a => n for (n, a) in enumerate(data_BH.tok_shots)])

	CH = [current_heating("aug", shot) for shot in pred_shots]
	
	col_dict = Dict("baseline"=>1, "hybrid"=>2, "training"=>3)
	CO_dict = Dict("CO"=> :cross, "C-EH"=> :ltriangle, "C-LH"=> :rtriangle)

	training_D = DataFrame(:ts => vcat(Bchoices, Hchoices),
						:predict => vcat(["training" for _ in 1:length(Bchoices)],
							["training" for _ in 1:length(Hchoices)]),
						:conf => vcat([(ts[2], (0.6*cos(n+1))) for (n, ts) in enumerate(Bchoices)],
									[(ts[2], (0.4*cos(n+1))) for (n, ts) in enumerate(Hchoices)]))
	training_D.CH = [current_heating(ts...) for ts in training_D.ts]

	RES.conf = [Tuple.(zip(ts[2], CONFIDENCE[1, n])) for (n, ts) in enumerate(RES.ts)]
	
	# RES_updated = vcat(RES[!, [:ts, :predict, :conf]], training_D)
	RES.CH .= CH
	ỹ = RES.predict
	
	HYBRID = (@subset RES @byrow :predict == "hybrid")[1:4:end, :]
	CO = (@subset RES @byrow :CH == "CO")
	
	f = Figure(size=(1200, 700));
	a0 = Axis(f[1:3, 1:12])
	a = Axis(f[4:16, 1:12])
	a1 = Axis(f[17:19, 1:12])
    a_sp = scatter!(a, RES.conf; colormap=:tab10, 
		color=[col_dict[el] for el in RES.predict], 
		colorrange = (1, 3),
		strokewidth = 0.1, 
		markersize = 18,
		marker = [CO_dict[label] for label in RES.CH],
		label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ])
	
		a0.xticklabelsvisible = false
	a0.yticklabelsvisible = false
	hidespines!(a0)
	a0.limits = (nothing, (-0.7, 0.7))
	hidedecorations!(a0)
	a.xticks=(10000:10000:40000, ["10,000", "20,000", "30,000", "40,000"])
	a1.yticklabelsvisible = false
	hidespines!(a1)
	a1.limits = (nothing, (-0.6, 0.6))
	hidedecorations!(a1)
	
	linkxaxes!(a0, a, a1)
	
	a0_sp = scatter!(a0, training_D.conf[1:length(Bchoices)]; colormap=:tab10, 
		color=3, 
		colorrange = (1, 3),
		strokewidth = 0.1, 
		markersize = 18,
		marker = [CO_dict[label] for label in training_D.CH[1:length(Bchoices)]],
		label = "training")

	scatter!(a1, training_D.conf[length(Bchoices)+1:end]; colormap=:tab10, 
		color=3, 
		colorrange = (1, 3),
		strokewidth = 0.1, 
		markersize = 18,
		marker = [CO_dict[label] for label in training_D.CH[length(Bchoices)+1:end]],
		label = "training"
	)

	pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
    for (n, (ts, coord)) in enumerate(zip(HYBRID.ts, HYBRID.conf))
        text!(a, coord, text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13)
    end
	RES_base_CO = @subset RES @byrow begin
		:predict == "baseline"
		:CH == "CO"
	end
    for (n, (ts, coord)) in enumerate(zip(RES_base_CO.ts, RES_base_CO.conf))
		text!(a, coord, text="$(ts[2])", align=pos[mod(n+1,4)+1], fontsize=13)
    end
	pos = [(:right, :bottom), (:left, :top), (:left, :bottom), (:right, :top)]
	tr_D_CO = @subset training_D @byrow :CH == "CO"
    for (n, (ts, coord)) in enumerate(zip(tr_D_CO.ts, tr_D_CO.conf))
		text!(a1, coord, text="$(ts[2])", align=pos[mod(n+2,4)+1], fontsize=13)
    end

	# extra = (@subset training_data @byrow in(:ts, [("aug", shot) for shot in [12000, 12032]]))
	# for (n, (ts, coord)) in enumerate(zip(extra.ts, extra.conf))
 #        text!(a, coord, text="$(ts[2])", align=pos[mod(n,2+2)+1], fontsize=7)
 #    end

	# a.xgridvisible = false
	# a.ygridvisible = false
	
	group_marker = [MarkerElement(marker = marker, color = :black,
	    strokecolor = :transparent,
	    markersize = 12) for marker in [:ltriangle, :rtriangle, :cross]]

	Legend(f[15, 2],
	    group_marker,
	    ["EH", "LH", "CO"],
	 	framevisible=false, labelsize=20)

	colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
	group_marker = [MarkerElement(marker = :circ, color = colors_3[i],
	    strokecolor = :transparent,
	    markersize = 12) for i in 1:3]

    Legend(f[15, 1],
	    group_marker,
	    ["baseline", "hybrid", "training"],
	 	framevisible=false, labelsize=20)
	# axislegend(a0, unique=true, merge=true, position=:lb, labelsize=20)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/Baseline_Hybrid_classified.png", f)
	
	# println("hybrid = ", (@subset RES @byrow :predict == "hybrid").CH |> countmap)
	# println("baseline = ", (@subset RES @byrow :predict == "baseline").CH |> countmap)
	
    f
end
let
	f = Figure(size=(1800, 900));
	a1 = Axis(f[1:3, 1:20], title="Baseline", ylabel="testing shots", ylabelsize=20, yticks=1:50:500, ylabelpadding=20)
	a2 = Axis(f[1:3, 21:25], yticklabelsvisible=false, yticksvisible=false, title="Hybrid")
	he = heatmap!(a1, KERNEL[1:length(Bchoices), :])
	heatmap!(a2, KERNEL[length(Bchoices)+1:end, :])
	Colorbar(f[1:3, end+1], he, ticks=0:0.1:1, width=20)
	a1.xticks = (1:length(Bchoices), ["#$j" for (i,j) in Bchoices])
	a2.xticks = (1:length(Hchoices), ["#$j" for (i,j) in Hchoices])
	
	a1.xticklabelrotation = π/2
	a2.xticklabelrotation = π/2
	
	a1.xticklabelsize = 15
	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/BH_kernel_matrix.png", f)
	f
end

begin
	# 36177, 29772, 13622, 34610, 17282, 15714, 19314, 32135, 34441, 32464, 33407
	shot = 19112
	x1_lim = 8.5
	x2_lim = 4
	fig1 = Figure(size=(1700, 800))
	ax1 = Axis(fig1[1:3, 2:4], title="#$shot", titlesize=25)
	ax2 = Axis(fig1[4:6, 2:4])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}"]))
		P = profiles(("aug", shot)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=:dracula, colorrange=(1, 8), color=n, label=label)
	end
	ax1.limits = ((0., x1_lim), (-0.1, 2.5))
	ax2.limits = ((0., x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)
	
	shot = 11126
	ax3 = Axis(fig1[1:3, 5:7], title="#$shot", titlesize=25)
	ax4 = Axis(fig1[4:6, 5:7])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1], [ax3, ax3, ax3, ax3, ax4, ax4, ax4], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}"]))
		P = profiles(("aug", shot)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=:dracula, colorrange=(1, 8), color=n, label=label)
	end
	ax3.limits = ((0., x2_lim), (-0.1, 2.5))
	ax4.limits = ((0., x2_lim), (-0.1, 2.5))
	ax3.yticks = 0:0.5:2.5
	ax4.yticks = 0:0.5:2.5
	# fig1[1:2, 8] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	# fig1[4:5, 8] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/shot_comparison.png", fig1)
	linkxaxes!(ax1, ax2)
	linkxaxes!(ax3, ax4)
	fig1
end

profiles("aug", 15714)