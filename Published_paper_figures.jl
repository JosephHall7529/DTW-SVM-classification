# Table 1
let
    AUG_data = (@subset all_2D_features @byrow begin
        :tok == "aug"
        :current_heating !== missing
    end)
    CM = cardinality_metadata(AUG_data, :current_heating)[1:3, :]
    TM = DataFrame(trial_dataset_meta)[2:end, :]
    insertcols!(TM, 1, :current_heating => ["C-LH", "C-EH", "CO"])
    trial_D = innerjoin(CM, TM, on=:current_heating)
    # trial_D |> clipboard
end

#Figure 1 
GLMakie.activate!()
begin
	signal_1 = [0, 0, 0, 0, 1, 1, 0, 0, 0, -2, -1, 0, 0]
	signal_2 = [0, 0, 1, 1, 0, 0, -2, 0, 0, 0, 0]
	mat = dtw_cost_matrix(signal_2, signal_1, Cityblock())
	cost, i1, i2 = DynamicAxisWarping.trackback(mat)
end
begin
    fig = Figure(size=(1000, 550));
    gl1 = fig[1, 1] = GridLayout()
    begin
        a11 = Axis(gl1[1, 1]; xticks=0:13, xlabel="time", yticklabelsvisible=false, yticksvisible=false, aspect = AxisAspect(2))
        hidedecorations!(a11, ticks=false, ticklabels=false)
        hidespines!(a11, :t, :l, :r)
        lines!(a11, signal_1 .+ 0.6, label="signal_1")
        lines!(a11, signal_2 .- 0.6, label="signal_2")
        
    end
    gl2 = fig[1, 2] = GridLayout()
    begin
        ax = Axis(gl2[1:11, 3:14])
        heatmap!(mat, colormap=:thermal)
        lines!(ax, i2, i1, color=:white, linewidth=3)
        begin
            ax.xticklabelsvisible=false
            ax.yticklabelsvisible=false
            ax.xticks=1:13
            ax.xticksvisible=false
            ax.yticksvisible=false
        end
    
        for (i, j) in Iterators.product(1:13, 1:11)
            text!((i-0.1, j-0.3),
                text="$(mat[i, j])",
                color=(:white, 0.5)
            )
        end
        ax1 = Axis(gl2[1:11, 1:2])
        hidedecorations!(ax1)
        hidespines!(ax1)
        lines!(ax1, -signal_2, 1:11, color=:orange)
    
        ax2 = Axis(gl2[12:14, 3:14])
        hidedecorations!(ax2)
        hidespines!(ax2)
        lines!(ax2, signal_1, color=:blue)
    end
    gl3 = fig[2, 1:2] = GridLayout()
    begin
        ax = Axis(gl3[1, 1], xticks=1:13, yticklabelsvisible=false, yticksvisible=false, aspect = DataAspect())
        hidespines!(ax, :t, :l, :r)
        hidedecorations!(ax, ticks=false, ticklabels=false)

        znorm(x) = (x = x.- mean(x); x ./= std(x))

        separation, ds = 0.7, 1
        x, y = (signal_2, signal_1)

        s1 = x .- separation
        s2 = y .+ separation

        i = fill(Inf, 1, length(i1))
        p1, p2 = vec([i1'; i2'; i][:, 1:ds:end]), vec([s1[i1]'; s2[i2]'; i][:, 1:ds:end])

        lines!(s2)
        lines!(s1)
        lines!(p1, p2, color=(:gray, 0.5))
    end
    # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/DTW_example.png", fig)
	fig
end

# Figure 2
let
	# 36177, 29772, 13622, 34610, 17282, 15714, 19314, 32135, 34441, 32464, 33407
    fig1 = Figure(size=(1000, 400))
	shot = 34841 
	
    x1_lim = 8.5
	x2_lim = 4
	
	ax1 = Axis(fig1[1:3, 2:4])
	ax2 = Axis(fig1[4:6, 2:4])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, color=param_colors[n], label=label)
	end
	ax1.limits = ((0., x1_lim), (-0.1, 2.5))
	ax2.limits = ((0., x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)
	
	shot = 34770

	ax3 = Axis(fig1[1:3, 5:7])
	ax4 = Axis(fig1[4:6, 5:7])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, color=param_colors[n], label=label)
	end
	ax3.limits = ((0., x2_lim), (-0.1, 2.5))
	ax4.limits = ((0., x2_lim), (-0.1, 2.5))
	ax3.yticks = 0:0.5:2.5
	ax4.yticks = 0:0.5:2.5
	# fig1[1:2, 8] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	# fig1[4:5, 8] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/shot_comparison.png", fig1)
	fig1
end

# Figure 3
let
	FTIP = flat_top_IP(tok_shots(which(["IP"])), 0.8)
    FTNBI = flat_top_NBI(tok_shots(which(["PNBI"])), 0.7)

	shot = 12021

	x1_lim = 6.6
	x2_lim = 4.5

	fig1 = Figure(size=(1000, 400))

	ax1 = Axis(fig1[1:3, 5:7])
	ax2 = Axis(fig1[4:6, 5:7])

	features = ["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"]
	t = FTIP[("aug", shot)]
    
	for (n, (feat, norm, axis, label)) in enumerate(zip(features, [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		BV = t[1] .< P.t .< t[2]
		Y = abs.(P.y.y .* norm)
		x_shade = P.t[BV]
		y_shade = Y[BV]
		band!(axis, x_shade, 0, y_shade, color=(:cyan, 0.15))
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, color=param_colors[n], label=label)
	end
	ax1.limits = ((0, x1_lim), (-0.1, 2.5))
	ax2.limits = ((0, x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)
	
	shot = 12021

	t = FTNBI[("aug", shot)]

	ax3 = Axis(fig1[1:3, 2:4])
	ax4 = Axis(fig1[4:6, 2:4])

	for (n, (feat, norm, axis, label)) in enumerate(zip(features, [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], 
                                                        [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], 
                                                        [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		BV = t[1] .< P.t .< t[2]
		Y = abs.(P.y.y .* norm)
		x_shade = P.t[BV]
		y_shade = Y[BV]
		band!(axis, x_shade, 0, y_shade, color=(:purple, 0.15))
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, color=param_colors[n], label=label)
	end
	ax3.limits = ((0, x1_lim), (-0.1, 2.5))
	ax4.limits = ((0, x1_lim), (-0.1, 2.5))
	ax3.yticks = 0:0.5:2.5
	ax4.yticks = 0:0.5:2.5
	# fig1[1:2, 8] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
	# fig1[4:5, 8] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)

	# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/12021_TIs.png", fig1)
	fig1
end

# Figure ?
let 
    CLn, CEn, COn = 20, 20, 10
    CLchoices = sort(sample(Random.seed!(13), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false))
    CEchoices = sort(sample(Random.seed!(13), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false))
    COchoices = sort(sample(Random.seed!(13), tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false))

    f = Figure();
    a = Axis(f[1, 1], xticks=0:5000:45000, xtickformat="{:.0f}", ygridvisible=false, yticklabelsvisible=false, yticksvisible=false)

    shots = vcat([el[2] for el in CLchoices],
        [el[2] for el in CEchoices], 
        [el[2] for el in COchoices]
    )
    markers = vcat([:ltriangle for _ in 1:20],
        [:rtriangle for _ in 1:20],
        [:cross for _ in 1:10]
    )

    y = vcat([0.1 for _ in 1:20],
        [0.2 for _ in 1:20],
        [0.3 for _ in 1:10]
    )

    scatter!(a, shots, y, marker=markers, markersize=15)

    f
end 

# Plot the train test case
let 
	# CairoMakie.activate!()
	f = Figure(size=(1000, 900))
	tests = [15714, 33371]
	opaque = ["IP", "PNBI"]
	σ = false
	begin 
        trial = 1
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
		train_ind = training_partion(labelled_data, unique(labelled_y), k, S=11)

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
        trial = 1
        para = naming(["BETAPOL"])
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
	train = 42432
	test = 16788
	
	x1_lim = 10
	x2_lim = 10
	
	compare_plot(("aug", train), ("aug", test), (0, 5), (0, 5); band="NBI")
end