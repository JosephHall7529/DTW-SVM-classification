function plot_layout(ℓ)
    @assert ℓ <= 10
    if ℓ <= 3
        layout = [(j, 1) for j in 1:ℓ]
        □ = (700, 500*ℓ)
    elseif ℓ == 5
        layout = [(j, 1) for j in 1:ℓ]
        □ = (700, 500*ℓ)
    elseif mod(ℓ, 2) == 0
        layout = [(i, j) for (i, j) in Iterators.product(1:div(ℓ, 2), 1:2)]
        □ = (1400, 500*div(ℓ, 2))
    else
        layout = vcat([(j, 1) for j in 1:div(ℓ, 2)+1], [(j, 2) for j in 1:div(ℓ, 2)])
        □ = (1400, 500*(div(ℓ, 2)+1))
    end
    return layout, □
end

# plotting variability of regression parameters
function regression_viability(D::DataFrame, dependent::String, predictors::Vector{String}; kwargs...)
    parameters = vcat(dependent, predictors)
    ℓ = length(parameters)
    layout, □ = plot_layout(ℓ)

    fig = Figure(size=□)

    axes = []
    for (k, sub_) in enumerate(layout)
        ax = Axis(fig[sub_...]; kwargs...)
        hist!(ax, D[:, parameters[k]], normalization=:probability)
        ax.xlabel = LaTeX_dict[parameters[k]]
        ax.ylabel = "P(X=x)"
        push!(axes, ax)
    end
    fig
end

# regression_viability(AUG, "TAUTH", ["IP", "BT", "NEL", "PLTH"])

# plotting time traces
function plot_1D_features(ts::Tuple{String, Number}, feat::String; kwargs...)
    P = profiles(ts..., feat)  

    fig, ax, line = lines(P.t, abs.(P.y[:, 1].*normalise_2D_features[feat]); kwargs...)
    ax.xlabel = "t"
    ax.ylabel = feat
    ax.title = "$ts"
    # axislegend(ax)
    return fig
end
function plot_1D_features(ts::Tuple{String, Number}, feats::Vector{String}; kwargs...)
    ℓ = length(feats)

    P = [profiles(ts..., feat) for feat in feats]  
    Y = [abs.(P[n].y.y.*normalise_2D_features[feat]) for (n, feat) in enumerate(feats)]

    fig = Figure();
    ax = Axis(fig[1, 1]);

    for (k, feat) in enumerate(feats)
        lines!(ax, P[k].t, Y[k], label=feat)
    end
    ax.xlabel = "t"
    axislegend(ax)
    ax.title = "$ts"

    return fig, ax
end
function plot_1D_features(ts::Tuple{String, Number}, feat::String, shade::Bool; short_IP::Bool=false, kwargs...) 
    if !shade
        return plot_1D_features(ts, feat)
    end 

    P = profiles(ts..., feat)
    Y = abs.(P.y.y .* normalise_2D_features[feat])

    if short_IP
        t = Dict("IP" => FTIPs[ts], "NBI" => FTNBI[ts])
    else
        t = Dict("IP" => FTIP[ts], "NBI" => FTNBI[ts])
    end 
    BV = Dict("IP" => (t["IP"][1] .< P.t .< t["IP"][2]), "NBI" => (t["NBI"][1] .< P.t .< t["NBI"][2]))

    x_shade = Dict("IP" => P.t[BV["IP"]], "NBI" => P.t[BV["NBI"]])
    y_shade = Dict("IP" => Y[BV["IP"]], "NBI" => Y[BV["NBI"]])

    fig = Figure();
    ax = Axis(fig[1, 1]; kwargs...)
    lines!(ax, P.t, Y)
    band!(ax, x_shade["IP"], 0, y_shade["IP"], color=(:red, 0.2), label="flat top IP")
    band!(ax, x_shade["NBI"], 0, y_shade["NBI"], color=(:blue, 0.2), label="flat top NBI")

    ax.xlabel = "t"
    ax.ylabel = feat
    ax.title = "$ts"
    Legend(fig[1, 2], ax, framevisible=false)
    return fig
end
function plot_1D_features(ts::Tuple{String, Number}, feats::Vector{String}, shade::Bool; short_IP::Bool=false, kwargs...)
    if !shade
        return plot_1D_features(ts, feats)
    end 
    
    ℓ = length(feats)
    layout, □ = plot_layout(ℓ) 

    P = [profiles(ts..., feat) for feat in feats]  
    Y = [abs.(P[n].y.y.*normalise_2D_features[feat]) for (n, feat) in enumerate(feats)]

    if short_IP
        t = Dict("IP" => FTIPs[ts], "NBI" => FTNBI[ts])
    else
        t = Dict("IP" => FTIP[ts], "NBI" => FTNBI[ts])
    end 
    BV = Dict("IP" => [(t["IP"][1] .< P[k].t .< t["IP"][2]) for k in 1:ℓ], "NBI" => [(t["NBI"][1] .< P[k].t .< t["NBI"][2]) for k in 1:ℓ])
    
    x_shade = Dict("IP" => [P[k].t[BV["IP"][k]] for k in 1:ℓ], "NBI" => [P[k].t[BV["NBI"][k]] for k in 1:ℓ])
    y_shade = Dict("IP" => [Y[k][BV["IP"][k]] for k in 1:ℓ], "NBI" => [Y[k][BV["NBI"][k]] for k in 1:ℓ])

    fig = Figure(size=□);
    
    axes = []
    for (k, sub_) in enumerate(layout)
        ax = Axis(fig[sub_...]; kwargs...)
        lines!(ax, P[k].t, Y[k])
        band!(ax, x_shade["IP"][k], 0, y_shade["IP"][k], color=(:red, 0.2), label="90% IP")
        band!(ax, x_shade["NBI"][k], 0, y_shade["NBI"][k], color=(:blue, 0.2), label="80% NBI")
        ax.xlabel = "t"
        ax.ylabel = "$(feats[k])"
        push!(axes, ax)
    end
	axes[1].title = "$(ts[1]) #$(ts[2])"
    axislegend(axes[minimum([4, ℓ])], position=:lt)
    display(fig)
    return axes, fig
end

function overview_plot(ts::Tuple{String, Number}, x1_lim::Union{Nothing, Tuple}=nothing; disp::Bool=false, band::String="")
    fig1 = Figure(size=(1000, 400))
	
	ax1 = Axis(fig1[1:3, 2:4], title="$(ts)")
	ax2 = Axis(fig1[4:6, 2:4])

	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], 
                                                            [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], 
                                                            [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], 
                                                            [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)",
                                                                L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(ts..., feat) 
        Y = abs.(P.y.y .* norm)
        if band != ""
            t = Dict("IP" => FTIP[ts], "NBI" => FTNBI[ts])
            BV = Dict("IP" => (t["IP"][1] .< P.t .< t["IP"][2]), "NBI" => (t["NBI"][1] .< P.t .< t["NBI"][2]))

            x_shade = Dict("IP" => P.t[BV["IP"]], "NBI" => P.t[BV["NBI"]])
            y_shade = Dict("IP" => Y[BV["IP"]], "NBI" => Y[BV["NBI"]])

            if band == "IP"
                band!(axis, x_shade["IP"], 0, y_shade["IP"], color=(:cyan, 0.2))
            elseif band == "NBI"
                band!(axis, x_shade["NBI"], 0, y_shade["NBI"], color=(:red, 0.2))
            end
        end
		lines!(axis, P.t, Y, linewidth=1.5, color=param_colors[n], label=label)
	end
    linkxaxes!(ax2, ax1)
	ax1.limits = (x1_lim, nothing)
	ax2.limits = (x1_lim, (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)
	
    # if disp
    #     # display(GLMakie.Screen(), fig1)
    # else
        fig1
    # end
end

function compare_plot(shot_1::Tuple{String, Int}, shot_2::Tuple{String, Int}, x1_lim::Union{Nothing, Tuple}=nothing, x2_lim::Union{Nothing, Tuple}=nothing; disp::Bool=false, band::String="")
    fig1 = Figure(size=(1500, 500))
	
    gl1 = fig1[1, 1] = GridLayout()
	begin
        ax1 = Axis(gl1[1:3, 2:4], title="$(shot_1)")
        ax2 = Axis(gl1[4:6, 2:4])

        for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], 
                                                                [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], 
                                                                [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], 
                                                                [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)",
                                                                    L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
            P = profiles(shot_1..., feat) 
            Y = abs.(P.y.y .* norm)
            if band != ""
                t = Dict("IP" => FTIP[shot_1], "NBI" => FTNBI[shot_1])
                BV = Dict("IP" => (t["IP"][1] .< P.t .< t["IP"][2]), "NBI" => (t["NBI"][1] .< P.t .< t["NBI"][2]))

                x_shade = Dict("IP" => P.t[BV["IP"]], "NBI" => P.t[BV["NBI"]])
                y_shade = Dict("IP" => Y[BV["IP"]], "NBI" => Y[BV["NBI"]])

                if band == "IP"
                    band!(axis, x_shade["IP"], 0, y_shade["IP"], color=(:cyan, 0.2))
                elseif band == "NBI"
                    band!(axis, x_shade["NBI"], 0, y_shade["NBI"], color=(:red, 0.2))
                end
            end
            lines!(axis, P.t, Y, linewidth=1.5, color=param_colors[n], label=label)
        end
        linkxaxes!(ax2, ax1)
        ax1.limits = (x1_lim, nothing)
        ax2.limits = (x1_lim, (-0.1, 2.5))
        ax1.yticks = 0:0.5:2.5
        ax2.yticks = 0:0.5:2.5
        gl1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
        gl1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)
    end    

    gl2 = fig1[1, 2] = GridLayout()
    begin
        ax3 = Axis(gl2[1:3, 1:3], title="$(shot_2)")
        ax4 = Axis(gl2[4:6, 1:3])
    
        for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], 
                                                                [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], 
                                                                [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], 
                                                                [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)",
                                                                    L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
            P = profiles(shot_2..., feat) 
            Y = abs.(P.y.y .* norm)
            if band != ""
                t = Dict("IP" => FTIP[shot_2], "NBI" => FTNBI[shot_2])
                BV = Dict("IP" => (t["IP"][1] .< P.t .< t["IP"][2]), "NBI" => (t["NBI"][1] .< P.t .< t["NBI"][2]))

                x_shade = Dict("IP" => P.t[BV["IP"]], "NBI" => P.t[BV["NBI"]])
                y_shade = Dict("IP" => Y[BV["IP"]], "NBI" => Y[BV["NBI"]])

                if band == "IP"
                    band!(axis, x_shade["IP"], 0, y_shade["IP"], color=(:cyan, 0.2))
                elseif band == "NBI"
                    band!(axis, x_shade["NBI"], 0, y_shade["NBI"], color=(:red, 0.2))
                end
            end
            lines!(axis, P.t, Y, linewidth=1.5, color=param_colors[n], label=label)
        end
        linkxaxes!(ax4, ax3)
        ax3.limits = (x2_lim, nothing)
        ax3.limits = (x2_lim, (-0.1, 2.5))
        ax4.yticks = 0:0.5:2.5
        ax4.yticks = 0:0.5:2.5
    end
    linkxaxes!(ax2, ax1)
    linkxaxes!(ax4, ax3)
    linkyaxes!(ax3, ax1)
    linkyaxes!(ax4, ax2)
    # if disp
    #     # display(GLMakie.Screen(), fig1)
    # else
    fig1
    # end
end
function compare_plot(shot_1::Tuple{String, Int}, shot_2::Tuple{String, Int}, x1_lim::Union{Nothing, Tuple}=nothing, x2_lim::Union{Nothing, Tuple}=nothing; disp::Bool=false, band::String="")
    fig1 = Figure(size=(1500, 500))
	
    gl1 = fig1[1, 1] = GridLayout()
	begin
        ax1 = Axis(gl1[1:3, 2:4], title="$(shot_1)")
        ax2 = Axis(gl1[4:6, 2:4])

        for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], 
                                                                [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], 
                                                                [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], 
                                                                [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)",
                                                                    L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
            P = profiles(shot_1..., feat) 
            Y = abs.(P.y.y .* norm)
            if band != ""
                t = Dict("IP" => FTIP[shot_1], "NBI" => FTNBI[shot_1])
                BV = Dict("IP" => (t["IP"][1] .< P.t .< t["IP"][2]), "NBI" => (t["NBI"][1] .< P.t .< t["NBI"][2]))

                x_shade = Dict("IP" => P.t[BV["IP"]], "NBI" => P.t[BV["NBI"]])
                y_shade = Dict("IP" => Y[BV["IP"]], "NBI" => Y[BV["NBI"]])

                if band == "IP"
                    band!(axis, x_shade["IP"], 0, y_shade["IP"], color=(:cyan, 0.2))
                elseif band == "NBI"
                    band!(axis, x_shade["NBI"], 0, y_shade["NBI"], color=(:red, 0.2))
                end
            end
            lines!(axis, P.t, Y, linewidth=1.5, color=param_colors[n], label=label)
        end
        linkxaxes!(ax2, ax1)
        ax1.limits = (x1_lim, nothing)
        ax2.limits = (x1_lim, (-0.1, 2.5))
        ax1.yticks = 0:0.5:2.5
        ax2.yticks = 0:0.5:2.5
        gl1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=15)
        gl1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=15)
    end    

    gl2 = fig1[1, 2] = GridLayout()
    begin
        ax3 = Axis(gl2[1:3, 1:3], title="$(shot_2)")
        ax4 = Axis(gl2[4:6, 1:3])
    
        for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], 
                                                                [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], 
                                                                [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], 
                                                                [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)",
                                                                    L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
            P = profiles(shot_2..., feat) 
            Y = abs.(P.y.y .* norm)
            if band != ""
                t = Dict("IP" => FTIP[shot_2], "NBI" => FTNBI[shot_2])
                BV = Dict("IP" => (t["IP"][1] .< P.t .< t["IP"][2]), "NBI" => (t["NBI"][1] .< P.t .< t["NBI"][2]))

                x_shade = Dict("IP" => P.t[BV["IP"]], "NBI" => P.t[BV["NBI"]])
                y_shade = Dict("IP" => Y[BV["IP"]], "NBI" => Y[BV["NBI"]])

                if band == "IP"
                    band!(axis, x_shade["IP"], 0, y_shade["IP"], color=(:cyan, 0.2))
                elseif band == "NBI"
                    band!(axis, x_shade["NBI"], 0, y_shade["NBI"], color=(:red, 0.2))
                end
            end
            lines!(axis, P.t, Y, linewidth=1.5, color=param_colors[n], label=label)
        end
        linkxaxes!(ax4, ax3)
        ax3.limits = (x2_lim, nothing)
        ax3.limits = (x2_lim, (-0.1, 2.5))
        ax4.yticks = 0:0.5:2.5
        ax4.yticks = 0:0.5:2.5
    end
    # linkxaxes!(ax2, ax1)
    # linkxaxes!(ax4, ax3)
    # linkyaxes!(ax3, ax1)
    # linkyaxes!(ax4, ax2)
    # if disp
    #     # display(GLMakie.Screen(), fig1)
    # else
    fig1
    # end
end
# begin
#     axes, fig = plot_1D_features(("aug", 20115), ["IP", "PNBI", "PECRH", "Q95", "BETAPOL"], false)
#     fig
# end
# using GLMakie
# GLMakie.activate!()
# axes, fig = plot_1D_features(("aug", 17220), features, limits=((0, 8), nothing))
# axes[4].limits = ((0, 8), (0, 1.3))
# # axes[3].limits = ((0, 4.0), nothing)

# plot_1D_features(("aug", 29185), ["IP", "PNBI"], true)


function DTW_plot(shot_1::Tuple{String, Int}, shot_2::Tuple{String, Int}, feat::String; section::String="IP")
    begin
        if section == "IP"
            dist = Cityblock()
            # dist = CosineDist()
            signal_1 = data_dtw_[feat].profile_data[shot_1]
            signal_2 = data_dtw_[feat].profile_data[shot_2]
        elseif section == "NBI"
            dist = Cityblock()
            signal_1 = data_dtw_[feat].flat_top_data[shot_1]
            signal_2 = data_dtw_[feat].flat_top_data[shot_2]
        end
        
        s1fn = size(signal_1, 1)
        s2fn = size(signal_2, 1)
        s1n = size(signal_1, 2)
        s2n = size(signal_2, 2)
        mat = dtw_cost_matrix(signal_2, signal_1, dist)
        cost, i1, i2 = DynamicAxisWarping.trackback(mat)
        println(cost)
    end
    begin
        fig = Figure(size=(1000, 550));
        gl1 = fig[1, 1] = GridLayout()
        begin
            x_vals = range(extrema(vcat(signal_1[s1fn, :], signal_2[s2fn, :]))..., length=11)
            a11 = Axis(gl1[1, 1]; xticks=round.(x_vals, digits=1), xlabel="time", yticklabelsvisible=false, yticksvisible=false, aspect = AxisAspect(2))
            hidedecorations!(a11, ticks=false, ticklabels=false)
            hidespines!(a11, :t, :l, :r)
            for i in 1:(s1fn-1)
                lines!(a11, signal_1[s1fn, :], signal_1[i, :], label="signal_1", color=colors[i])
                lines!(a11, signal_2[s2fn, :], signal_2[i, :], label="signal_2", color=colors[i])
            end
            
        end
        gl2 = fig[1, 2] = GridLayout()
        begin
            ax = Axis(gl2[1:11, 3:14])
            hm = heatmap!(mat, colormap=:thermal, colorscale=log10)
            Colorbar(gl2[1:11, 15:16], hm, scale=log10)
            lines!(ax, i2, i1, color=:white, linewidth=3)
            begin
                ax.xticklabelsvisible=false
                ax.yticklabelsvisible=false
                # ax.xticks=1:13
                ax.xticksvisible=false
                ax.yticksvisible=false
            end
        
            # for (i, j) in Iterators.product(1:s1n, 1:s2n)
            #     text!((i-0.1, j-0.3),
            #         text="$(mat[i, j])",
            #         color=(:white, 0.5)
            #     )
            # end
            ax1 = Axis(gl2[1:11, 1:2])
            hidedecorations!(ax1)
            hidespines!(ax1)
            for i in 1:s2fn-1
                lines!(ax1, -signal_2[i, :], 1:s2n, color=colors[i])
            end
        
            ax2 = Axis(gl2[12:14, 3:14])

            hidedecorations!(ax2)
            hidespines!(ax2)
            for i in 1:s1fn-1
                lines!(ax2, 1:s1n, signal_1[i, :], color=colors[i])
            end
        end
        gl3 = fig[2, 1:2] = GridLayout()
        begin
            ax = Axis(gl3[1, 1], xticks=0:10:s1n, yticklabelsvisible=false, yticksvisible=false)
            hidespines!(ax, :t, :l, :r)
            hidedecorations!(ax, ticks=false, ticklabels=false)

            znorm(x) = (x = x.- mean(x); x ./= std(x))

            separation, ds = -0.2, 5
            for j in 1:s1fn-1
                x, y = (signal_2[j, :], signal_1[j, :])

                s1 = x .- separation
                s2 = y .+ separation

                i = fill(Inf, 1, length(i1))
                p1, p2 = vec([i1'; i2'; i][:, 1:ds:end]), vec([s1[i1]'; s2[i2]'; i][:, 1:ds:end])

                lines!(s2, color=colors[j])
                lines!(s1, color=colors[j])
                lines!(p1, p2, color=(:gray, 0.5))
            end
        end
        # save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/DTW_example.png", fig)
    end
    display(GLMakie.Screen(), fig)
end

function overview_fit(model_::DTW_SVM, svmmodel, test_X, test_y=[""])
	f = Figure(size=(1000, 900))
	begin 
        fig_model = deepcopy(model_)
        training = fig_model.database.training_data

        labels = String.(unique(training.metadata.y))
        label_length = length(labels)

		MLJBase.predict(fig_model, svmmodel, test_X)
        testing = fig_model.database.testing_data

        if test_y !== [""]
            testing.metadata.actual = [test_y[ind] for ind in testing.metadata.order]
            sort!(testing.metadata, [:actual, :tok_shots], rev=[true, false])
        end
    end
    # explicitly calculating the kernel element of cosine and flat top
    begin
        df_names = df_ts_naming(testing.metadata.tok_shots)

		dat = Array(testing.cosine_cost[:, df_names])
		mat_c = exp.(-(dat ./ fig_model.C_cosine).^2)

		dat = Array(testing.flat_top_cost[:, df_names])
		mat_ft = exp.(-(dat ./ fig_model.C_flat_top).^2)

		mat = mat_c .* mat_ft
	end

	gl1 = f[1:2, 1] = GridLayout()
	begin
        col_dict = Dict("LH"=>3, "EH"=>2, "CO"=>1)

		axmain = Axis(gl1[1, 1], xlabel="training", ylabel="test")
		he = heatmap!(axmain, mat, colormap=:oslo)

		if test_y !==[""]
            cntmap = Dict([label => count(i -> i == label, testing.metadata.actual) for label in labels])
        else
            cntmap = Dict([label => count(i -> i == label, testing.metadata.y) for label in labels]) 
        end
        k = Dict([label => count(i -> i == label, training.metadata.y) for label in labels])
        labels = vcat([""], labels)
        cntmap[""] = 0
        k[""] = 0

        for i in 3:label_length+1
            vertical_x0 = [sum(k[labels[j]] for j in 1:i-1)+0.5]
            vertical_y0 = sum(cntmap[labels[j]] for j in 1:i-2) + 0.5 
            vertical_y1 = sum(cntmap[labels[j]] for j in 1:i) + 0.5

            horizontal_x0 = sum(k[labels[j]] for j in 1:i-2) + 0.5
            horizontal_x1 = sum(k[labels[j]] for j in 1:i) + 0.5
            horizontal_y0 = [sum(cntmap[labels[j]] for j in 1:i-1)+0.5] 
            
            lines!(axmain, vertical_x0, [vertical_y0, vertical_y1], color=:white, linewidth=4)
            lines!(axmain, [horizontal_x0, horizontal_x1], horizontal_y0, color=:white, linewidth=4)
        end

        labels = labels[2:end]
        
		Colorbar(gl1[1, 2], he, ticks=0:0.1:1)
	
		axmain.xticks = (1:length(training.metadata.tok_shots), ["#$j" for (i, j) in training.metadata.tok_shots])
		axmain.xticklabelrotation = π/2
		axmain.yticks = (1:length(testing.metadata.tok_shots), ["#$j" for (i, j) in testing.metadata.tok_shots])

		# axmain.title = "C-LH 								C-EH 									CO"
	end
	
	gl2 = f[1:3, 2] = GridLayout()
	let
        start = 1
        stop = 0
        for (n, label) in enumerate(labels)
            axleftCL = Axis(gl2[n, 1], xticks=(0.7:0.2:0.7, [L"\mathrm{DTW_{cos}}|_{I_p^{80%}}"]), xticklabelsize=20, 
                limits=(nothing, (0, 1.)), ylabel=label, ylabelsize=15)
            axmidCL = Axis(gl2[n, 2], xticks=(0.7:1:0.7, [L"\mathrm{DTW_{mag}}|_{P_{NBI}^{70%}}"]), xticklabelsize=20,
                limits=(nothing, (0, 1.)))
            axrightCL = Axis(gl2[n, 3], xticks=(0.7:1:0.7, [L"K"]), xticklabelsize=20, limits=(nothing, (0, 1.)))

            ℓ = size(training.metadata, 1)

            stop += k[label]
            rng = start:stop
            
            CD = vcat([mat_c[i, 1:end] for i in rng]...)
            hist!(axleftCL, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
            CD = vcat([mat_ft[i, 1:end] for i in rng]...)
            hist!(axmidCL, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
            CD = vcat([mat[i, 1:end] for i in rng]...)
            hist!(axrightCL, CD, normalization=:probability, scale_to=-0.6, offset=1, direction=:x, color=:gray70)
		end
        # if σ
		# 	for (N, te) in enumerate(teind) 
		# 		val, ind = findmax(mat[rng, te])
		# 		println("C-LH train_ind: ", labelled_ts[train_ind][rng[1]+ind-1], ", K = $(round(val, digits=2))")
		# 		hlines!(axleftCL, mat_c[rng[1]+ind-1, te], color=colors[N])
		# 		hlines!(axmidCL, mat_ft[rng[1]+ind-1, te], color=colors[N])
		# 		hlines!(axrightCL, val, color=colors[N])
		# 	end
		# end
	end
	
	gl3 = f[3, 1] = GridLayout()
	let
        ỹ = testing.metadata.y

        if test_y !== [""]
            ŷ = testing.metadata.actual

            println("Accuracy: ", Accuracy()(ỹ, ŷ), 
                "\nBalaced Accuracy: ", BalancedAccuracy()(ỹ, ŷ), 
                "\nFalse Positive Rate: ", MulticlassFalsePositiveRate(;levels=["LH", "EH", "CO"], perm=[3,2,1])(ỹ, ŷ))
        end
        conf = hcat([reshape(conf, :, 1) for conf in testing.metadata.confidence]...)

		axmain = Axis3(gl3[1:2, 1:5], azimuth=-0.53π, xlabelvisible=false, ylabelvisible=false,
			zlabelvisible=false, viewmode=:stretch, elevation=0.05π)
		axleg = Axis(gl3[1:2, 0])
		hidedecorations!(axleg)
		hidespines!(axleg)

        if test_y !== [""]
            scatter!(axmain, conf[1, :], conf[2, :], conf[3, :]; colormap=:tab10,
                color=[col_dict[el] for el in ỹ],
                strokewidth = 0.1,
                label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
                markersize = [(i == j ? 20 : 10) for (i, j) in zip(ŷ, ỹ)],
                marker = [(i == j ? '∘' : (:utriangle)) for (i, j) in zip(ŷ, ỹ)]
            )	
        else
            scatter!(axmain, conf[:, 1], conf[:, 2], conf[:, 3]; colormap=:tab10,
                color=[col_dict[el] for el in ỹ],
                strokewidth = 0.1,
                label = [label => (;colormap=:tab10, colorrange=(1, 3), color=col_dict[label]) for label in ỹ],
                markersize = 10,
                marker = utriangle
            )	
        end
		scatter!(axmain, [0], [0], [0], color=(:gray40, 0.3), markersize=20, strokewidth=0.2)
		
		pos = [(:right, :bottom), (:left, :bottom), (:right, :top), (:left, :top)]
		# shuffle!(Random.seed!(3), testing.metadata)
        for (n,row) in enumerate(eachrow(testing.metadata))
            text!(axmain, row.confidence..., text="$(row.tok_shots[2])", align=pos[mod(n+1,4)+1], fontsize=12, color=:gray40)
        end

		group_marker = [MarkerElement(marker = marker, color = :black,
			strokecolor = :transparent, markersize = markersize) for (marker, markersize) in zip([:utriangle, '∘'], [12, 27])]

		Legend(gl3[1, 0], group_marker, ["incorrectly \n predicted", "correctly \n predicted"], 
			framevisible=false, labelsize=14, halign=:left, valign=:bottom)

		colors_3 = get(ColorSchemes.tab10, range(0.0, 1.0, length=3))
		group_marker = [MarkerElement(marker = :circle, color = colors_3[i],
			strokecolor = :transparent, markersize = 12) for i in 1:3]

		Legend(gl3[2, 0], group_marker, ["LH", "EH", "CO"],
			framevisible=false, labelsize=14, halign=:left, valign=:top)

		axmain.xgridvisible = false
		axmain.ygridvisible = false
		axmain.zgridvisible = false
	end
	# colsize!(gl3, 1, Relative(0.05))
	
	# # Label(f[0, :], text = "", fontsize = 20)
	# # # save("/Users/joe/Project/PhD/EuroFusion/EUROfusion_ML_2024/Presentations/Multidimensional_qunatity_classification_20_01_25/CO_EH_LH_classification.png", f)
	# # display(GLMakie.Screen(), f)
	f
end

# function evaluation_plot(data_::Dict{Symbol, Dict{String, classification_performance}})
#     f = Figure();
#     a = Axis(f[1,1])
    
#     for i in 1:5
#         hist!(ax, randn(1000), scale_to=-0.6, offset=i, direction=:x)
#     end
# end