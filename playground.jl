select(AUG, [:id, :TOK, :SHOT, :TIME, :TAUTH, :τn98, :fbs, :Q95, :IP, :BT, :NEL, :nnorm, :RHOSTAR, :NUSTAR, :BETASTAR], :) |> vscodedisplay 

acc, K, cos_mid, mag_mid, res = hyper_parameter_search(data)


B = [i*j for i in 1:5, j in 1:5]
# upp_diag =
vcat([[B[i, j] for (i, j) in Iterators.product(k, k+1:5)] for k in 1:4]...)

res1 = (@subset res @byrow :label != "UNKNOWN")
confusion_matrix(res1.label, res1.predict)

profiles(("jet", 87213)...)
profiles(("aug", 20115)...)

JET = unique((@subset original_space(DATA_) @byrow in(:TOK, ["JET", "JETILW"])), :SHOT) |> vscodedisplay
JET_data = (@subset original_space(DATA_) @byrow begin
    in(:TOK, ["JET", "JETILW"])
    # :HYBRID == "YES"
    83630 .< :SHOT .< 92471
end)
select(JET_data, [:id, :TOK, :SHOT, :TIME, :TAUTH, :τn98, :fbs, :Q95, :IP, :BT, :NEL, :nnorm, :RHOSTAR, :NUSTAR, :BETASTAR], :) |> vscodedisplay 

AUG = unique((@subset original_space(DATA_) @byrow in(:TOK, ["AUG", "AUGW"])), :SHOT)
Regression(AUG, ols(), single_machine())

fs = split.(files("/Users/joe/Downloads/transfer_421294_files_e2609c28/jetdata_2025-03-05/data", false), "_")
fss = [emp.(String(el[1][2:end])) for el in fs] |> unique 
rest = setdiff(fss, JET_shots)  |> sort


# p = CSV.read("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/play.csv", DataFrame, stringtype=String)
for i in 1:size(res, 1)
    ts = res[i, :ts]
    int = @view p[(p.tok .== ts[1]) .& (p.shot .== ts[2]), :]
    int.current_heating .= res[i, :predict]
end
p
# CSV.write("/Users/joe/Project/Coding_clean/J_Hybrid_plasma_classification_25_02_24/data/play.csv", p)
.






shot_dict = Dict([a => n for (n, a) in enumerate(data.tok_shots)])
cos_step, mag_step = 5, 5
cos_mid, mag_mid = 50.5, 300.5


# CO -- EH
# C -- EH
# C -- LH
training_data = [27930, 34774, 42432, 26338,
    23008, 26519, 16746,
    16785, 29772, 27980,
    39124, 17151]
train_ts = [("aug", shot) for shot in training_data]
y = ["CO","CO","CO","CO",
    "C-EH","C-EH","C-EH",
    "C-LH","C-LH","C-LH",
    "C-PP", "C-PP"]

hyp_search, δ_cos, count, curr_acc = zeros(2, 2, 2), 1, 1, 0

while (count < 13) && (curr_acc !== 1.0) && (δ_cos > 0.01)
    Consts = vcat(collect(range(maximum([0.01, cos_mid-10*cos_step]), cos_mid, length=10)), [cos_mid + i*cos_step for i in 1:10])
    m = length(Consts)
    hyp_search = zeros(m, 2)

    for (N, CC) in enumerate(collect(Consts))

        X = exp.(-Array(data.cosine_cost[!, 2:end] ./ CC) .^2) .* exp.(-Array(data.magnitude_cost[!, 2:end] ./ 217).^2)

        TD_ind = vcat(sample(1:4, 1, replace=false),
            sample(5:7, 1, replace=false),
            sample(8:10, 1, replace=false),
            sample(11:12, 1, replace=false))
        TD = [(in(i, TD_ind) ? true : false) for i in 1:length(training_data)]
    
        ind = [shot_dict[k] for k in train_ts]
        
        K = X[ind[TD], ind[TD]]
        
        model = svmtrain(K, y[TD], kernel=LIBSVM.Kernel.Precomputed)
        
        KK = X[ind[TD], ind[Not(TD)]]
        ỹ, _ = svmpredict(model, KK)
    
        ACC = sum(y[Not(TD)] .== ỹ) / length(ỹ)
        hyp_search[N, 1] = ACC
        hyp_search[N, 2] = CC
    end
    max = maximum(hyp_search[:, 1])
    ind = findfirst(i -> i == max, hyp_search[:, 1])

    last_cos_mid = cos_mid
    (curr_acc, cos_mid) = (hyp_search[ind, 1], hyp_search[ind, 2])
    println(max, ": ($cos_mid)")
    δ_cos = abs(-(last_cos_mid, cos_mid))

    cos_step /= count+1

    count += 1
end
X = exp.(-Array(data.cosine_cost[!, 2:end] ./ cos_mid) .^2) .* exp.(-Array(data.magnitude_cost[!, 2:end] ./ 217).^2)
ind = [shot_dict[k] for k in train_ts]
K = X[ind, ind]

model = svmtrain(K, y, kernel=LIBSVM.Kernel.Precomputed)
KK = X[ind, Not(ind)]
ỹ, _ = svmpredict(model, KK)

res = DataFrame(hcat(data.tok_shots[Not(ind)], ỹ), [:ts, :predict])







# CO, C-EH ? 32103, 32104

shots = [j for (i, j) in res.ts]
n=1
n -= 1
CairoMakie.activate!()
begin
    shot = shots[n]
	println(shot)
	P = profiles(("aug", shot)..., "IP")
	fig, ax, lin = lines(P.t, P.y.y.*1e-6)
	P = profiles(("aug", shot)..., "PNBI")
	lines!(P.t, P.y.y.*1e-6)
	
	println(n)
    n+=1
	
	fig
end



BH_features = ["IP", "PNBI", "PECRH", "BETAPOL", "Q95", "NGW"]
DTW_hyp_1(tss, BH_features, 10, L=100)



18*0.6

data_try = let
    CEL_features = ["PICRH", "LI"]
    CEL_tss = tok_shots((@subset which(CEL_features) @byrow in(:tok, ["aug"])))
    DTW_hyp_1(CEL_tss, CEL_features, 10, L=100)
end

labelled_dict = let
    labelled_shots = vcat(Bchoices, Hchoices)
    labels = vcat(repeat(["baseline"], length(Bchoices)),
        repeat(["hybrid"], length(Hchoices))) 
    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
end

labelled_ts = collect(labelled_dict |> keys)
labelled_y = collect(labelled_dict |> values)
shot_dict = Dict([a => n for (n, a) in enumerate(data_try.tok_shots)])
labelled_ind = [shot_dict[k] for k in labelled_ts]

train_ind = training_partion(labelled_dict, ["baseline", "hybrid"], [37, 10])
CC = 6.
FTC = 160.

X = exp.(-(Array(data_try.cosine_cost[!, 2:end]) ./ CC).^2) .* 
    exp.(-(Array(data_try.flat_top_cost[!, 2:end]) ./ FTC).^2)

K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed)

KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
ỹ, conf = svmpredict(model, KK)
BalancedAccuracy()(ỹ, labelled_y[Not(train_ind)])

let
    col_dict = Dict("baseline"=>1, "hybrid"=>2, "training"=>3)
    pred_shots = [j for (i, j) in [data_try.tok_shots[i] for i in labelled_ind[Not(train_ind)]]]
    scatter(pred_shots, conf[1, :]; colormap=:tab10, 
        color=[col_dict[el] for el in ỹ], 
        colorrange = (1, 3),
        strokewidth = 0.1, 
        markersize = 18,
        marker = [pred == act ? '✓' : :cross for (pred, act) in zip(ỹ, labelled_y[Not(train_ind)])])
end




# ỹ, _ = KNN(KK, labelled_y[train_ind])





        
Regression()









begin
	function training_partion(labelled_data::OrderedDict, labels::Vector{String}, k::Int=2)
		ℓ = countmap(labelled_data |> values)

		i = 0
		TD_ind = Vector{Int}(undef, 0)
		for label in labels
			TD_ind = vcat(TD_ind, sample(i+1:i+ℓ[label], k, replace=false))
			i += ℓ[label]
		end
		
	    training_ind = [(in(i, TD_ind) ? true : false) for i in 1:length(labelled_data)]
	end
end
begin
    Random.seed!(17)
    CLn = 20
    CEn = 20
    COn = 5
    CLchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-LH"]))), CLn, replace=false)
    CEchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["C-EH"]))), CEn, replace=false)
    COchoices = sample(tok_shots((@subset all_2D_features @byrow in(:current_heating, ["CO"]))), COn, replace=false)
end
res, confidence, kernel_matrix = let
	labelled_data = let
	    labelled_shots = vcat(CLchoices, CEchoices, COchoices)
	    labels = vcat(repeat(["C-LH"], length(CLchoices)),
	        repeat(["C-EH"], length(CEchoices)),
	        repeat(["CO"], length(COchoices))) 
	    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
	end
	hyper_parameter_search(data, labelled_data, [4,4,2])
end
shots = [j for (i, j) in res.ts]
labels = (@subset all_2D_features @byrow begin
    :tok == "aug"
    in(:shot, shots)
end).current_heating
res.label .= labels

Accuracy()(labels, res.predict)
BalancedAccuracy()(labels, res.predict)








using GLMakie
GLMakie.activate!()
begin
	# 36177, 29772, 13622, 34610, 17282, 15714, 19314, 32135, 34441, 32464, 33407
	shot = 34841 
	x1_lim = 8.5
	x2_lim = 4
	fig1 = Figure(size=(1900, 900))
	ax1 = Axis(fig1[1:3, 2:4], title="#$shot", titlesize=25)
	ax2 = Axis(fig1[4:6, 2:4])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax1, ax1, ax1, ax1, ax2, ax2, ax2, ax2], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
		P = profiles(("aug", shot)..., feat)
		lines!(axis, P.t, abs.(P.y.y .* norm), linewidth=1.5, colormap=:dracula, colorrange=(1, 8), color=n, label=label)
	end
	ax1.limits = ((0., x1_lim), (-0.1, 2.5))
	ax2.limits = ((0., x1_lim), (-0.1, 2.5))
	ax1.yticks = 0:0.5:2.5
	ax2.yticks = 0:0.5:2.5
	fig1[1:2, 1] = Legend(fig1, ax1, "", framevisible=false, labelsize=35)
	fig1[4:5, 1] = Legend(fig1, ax2, "", framevisible=false, labelsize=35)
	
	shot = 34770
	ax3 = Axis(fig1[1:3, 5:7], title="#$shot", titlesize=25)
	ax4 = Axis(fig1[4:6, 5:7])
	for (n, (feat, norm, axis, label)) in enumerate(zip(["IP", "PNBI", "PECRH", "PICRH", "Q95", "BETAPOL", "NGW", "LI"], [1e-6, 1e-7, 1e-7, 1e-7, 1e-1, 1, 1, 1], [ax3, ax3, ax3, ax3, ax4, ax4, ax4, ax4], [L"I_p \, (10^{-6} MA)", L"P_{\mathrm{NBI}} \, (10^{-7} MA)", L"P_{\mathrm{ECRH}} \, (10^{-7} MA)", L"P_{\mathrm{ICRH}} \, (10^{-7} MA)", L"q_{95} \, (10^{-1})", L"\beta_p", L"f_{\mathrm{GW}}", L"\ell_i"]))
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
	fig1
end
# save("/Users/joe/Project/Papers_clean/Hybrid_classification_TH_T_24_02_22/figures/shot_comparison_stochastic.png", fig1)
let
	train_tss = vcat(CLchoices, CEchoices, COchoices)
	n = length(train_tss)
	ind = [data.shot_dict[ts] for ts in train_tss]

	fig = Figure();
	ax = Axis(fig[1, 1], yscale=sqrt)

    dat = data.cosine_cost
    mat = Array(dat[ind, ind.+1])
    CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
    hist!(ax, CD, scale_to=-0.6, offset=0, direction=:x)

    # dat = data.magnitude_cost
    # mat = Array(dat[ind, ind.+1])
    # CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
    # # println(describe(CD))
    # hist!(ax, CD, scale_to=-0.6, offset=1, direction=:x)

    # dat = data.flat_top_cost
    # mat = Array(dat[ind, ind.+1])
    # CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
    # # println(describe(CD))
    # hist!(ax, CD, scale_to=-0.6, offset=2, direction=:x)

	fig
end
let
	train_tss = vcat(CLchoices, CEchoices, COchoices)
	n = length(train_tss)
	ind = [data.shot_dict[ts] for ts in train_tss]

	fig = Figure();
	ax = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])

    dat = data.cosine_cost
    mat = Array(dat[ind, ind.+1])
    CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
    hist!(ax, CD, scale_to=-0.6, offset=0, direction=:x)

    dat = data.magnitude_cost
    mat = Array(dat[ind, ind.+1])
    CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
    # println(describe(CD))
    hist!(ax2, CD, scale_to=-0.6, offset=1, direction=:x)

    dat = data.flat_top_cost
    mat = Array(dat[ind, ind.+1])
    CD = vcat([mat[i, i+1:n] for i in 1:n-1]...)
    # println(describe(CD))
    hist!(ax2, CD, scale_to=-0.6, offset=2, direction=:x)

	fig
end
























Bchoices
Hchoices

BH_features = ["IP", "PNBI", "PICRH", "PECRH", "BETAPOL", "Q95", "NGW", "LI"]
BH_tss = tok_shots((@subset which(BH_features) @byrow in(:tok, ["aug"])))
data_main = DTW_hyp_1(BH_tss, BH_features, 10, L=100)

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
data_main
labelled_data = let
    labelled_shots = vcat(Bchoices, Hchoices)
    labels = vcat(repeat(["baseline"], length(Bchoices)),
        repeat(["hybrid"], length(Hchoices))) 
    d = OrderedDict([ts for ts in labelled_shots] .=> labels)
end
k = [15, 9]
interesting = "hybrid"
N = 10

labelled_ts = collect(labelled_data |> keys)
labelled_y = collect(labelled_data |> values)
shot_dict = Dict([a => n for (n, a) in enumerate(data_main.tok_shots)])
labelled_ind = [shot_dict[k] for k in labelled_ts]

# mag_arr = Array(data.magnitude_cost[:, 2:end])
ft_arr = Array(data_main.flat_top_cost[:, 2:end])
for i in 1:size(ft_arr, 2)
    ft_arr[i, i] = mean(ft_arr)
end
ft_range = extrema(ft_arr)

cos_arr = Array(data_main.cosine_cost[!, 2:end])
for i in 1:size(cos_arr, 2)
    cos_arr[i, i] = mean(cos_arr)
end
cos_range = extrema(cos_arr)

# cos_med = median(cos_arr[labelled_ind, labelled_ind][:])
# mag_med = median(mag_arr[labelled_ind, labelled_ind][:])
# ft_med = median(ft_arr[labelled_ind, labelled_ind][:])

shot_dict = Dict([a => n for (n, a) in enumerate(data_main.tok_shots)])
labelled_ts = collect(labelled_data |> keys)
labelled_y = collect(labelled_data |> values)
labels = collect(labelled_data |> values) |> unique

cnt, ACC, δ_acc, no_good_values = 1, 0., 100, 100

# while (cnt < 12) && (ACC !== 1.) && (δ_acc > 0.001)

# xs = range(ft_range..., length=3)
# ys = range(cos_range..., length=5)
xs = range(5, 10, length=100)
ys = range(77.5, 80, length=100)
Cs = collect(Iterators.product(
    xs,
    # range(mag_med/2, mag_med*2, length=10),
    ys
))

nc, nft = size(Cs)
hyp_search = zeros(nc, nft, 3)
acc_int = [zeros(nc, nft) for _ in 1:nthreads()]

Threads.@threads for S in 1:N
    Random.seed!(S)
    train_ind = training_partion(labelled_data, labels, k)
    for (i, j) in Iterators.product(1:nc, 1:nft)
        CC = Cs[i, j][1]
        FTC = Cs[i, j][2]

        X = exp.(-(Array(data_main.cosine_cost[!, 2:end]) ./ CC).^2) .* 
            exp.(-(Array(data_main.flat_top_cost[!, 2:end]) ./ FTC).^2) 

        K = X[labelled_ind[train_ind], labelled_ind[train_ind]]
        model = svmtrain(K, labelled_y[train_ind], kernel=Kernel.Precomputed)

        KK = X[labelled_ind[train_ind], labelled_ind[Not(train_ind)]]
        ỹ, _ = svmpredict(model, KK)
        # ỹ, _ = KNN(KK, labelled_y[train_ind])
        acc_int[threadid()][i, j] += BalancedAccuracy()(labelled_y[Not(train_ind)], ỹ)
        hyp_search[i, j, 2] = CC
        hyp_search[i, j, 3] = FTC
    end
end
hyp_search[:, :, 1] = map(+, acc_int...)
hyp_search[:, :, 1] ./= N
# hyp_search |> display

using CairoMakie
fig, ax, he = heatmap(xs, ys, hyp_search[:, :, 1]);
ax.xlabel = "flat top"
ax.ylabel = "cosine dist"
Colorbar(fig[:, end+1], he)

for (j, i) in Iterators.product(1:11:length(ys), 1:11:length(xs))
    text!((xs[i], ys[j]),
        text="$(round(hyp_search[i, j, 1], digits=4))",
        # text = "$(Cs[i, j])",
        color=(:gray, 0.5)
    )
end

fig



    max = maximum(hyp_search[:, :, 1])
    inds = findall(i -> i == max, hyp_search[:, :, 1])
    ind = sample(inds, 1)[1]

    # last_cos_med, last_mag_med = cos_med, mag_med
    (cos_med, ft_med) = (hyp_search[ind, 2], hyp_search[ind, 3])
    println(max, ": ($cos_med, $ft_med)")
    δ_acc = abs(-(ACC, max))

    ACC = max    
    # no_good_values = length(inds)
    cnt += 1
# end
X = exp.(-(Array(data.cosine_cost[!, 2:end]) ./ (cos_med)).^2) .* 
    exp.(-(Array(data.flat_top_cost[!, 2:end]) ./ (ft_med)).^2) 

K = X[labelled_ind, labelled_ind]
model = svmtrain(K, labelled_y, kernel=Kernel.Precomputed)

KK = X[labelled_ind, Not(labelled_ind)]
# ỹ, _ = KNN(KK, labelled_y)
ỹ, confidence = svmpredict(model, KK)

res = DataFrame(hcat(
        data.tok_shots[Not(labelled_ind)], ỹ
    ), [:ts, :predict]
)
if interesting !== ""
    show(stdout, "text/plain", (@subset res @byrow :predict == interesting))
end


