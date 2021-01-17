## Set up environment and load required packages
using NGWP, JLD, MAT, Plots, LightGraphs, Distances, MTSG
gr(dpi = 200)

## Build weighted toronto street network graph
G = loadgraph(joinpath(@__DIR__, "..", "datasets", "new_toronto_graph.lgz")); N = nv(G)
X = load(joinpath(@__DIR__, "..", "datasets", "new_toronto.jld"),"xy")
W = 1.0 .* adjacency_matrix(G)
dist_X = pairwise(Euclidean(),X; dims = 1)
Weight = W .* dualGraph(dist_X; method = "inverse") # weighted adjacence matrix
L = Matrix(Diagonal(sum(Weight;dims = 1)[:]) - Weight)
𝛌, 𝚽 = eigen(L); sgn = (maximum(𝚽, dims = 1)[:] .> -minimum(𝚽, dims = 1)[:]) .* 2 .- 1; 𝚽 = 𝚽 .* sgn';
Q = incidence_matrix(G; oriented = true)
edge_weight = 1 ./ sqrt.(sum((Q' * X).^2, dims = 2)[:])

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝚽,Q,N; edge_weight = edge_weight)
W_dual = sparse(dualGraph(distDAG)) #required: sparse dual weighted adjacence matrix

## Assemble natural graph wavelet packets
ht_elist_PC, ht_vlist_PC = HTree_EVlist(𝚽,W_dual)
wavelet_packet_PC = HTree_wavelet_packet(𝚽,ht_vlist_PC,ht_elist_PC)
ht_elist_VM = ht_elist_PC
wavelet_packet_VM = HTree_wavelet_packet_varimax(𝚽,ht_elist_VM)

## Fig. 10(a) a smooth spatial distribution of the street intersections graph signal
f = zeros(N); for i in 1:N; f[i] = length(findall(dist_X[:,i] .< 1/minimum(edge_weight))); end #fneighbor
gplot(W, X; width=1); scatter_gplot!(X; marker = f, smallValFirst = true, ms = 3); signal_plt = plot!(aspect_ratio = 1, grid = false)
savefig(signal_plt, joinpath(@__DIR__, "../paperfigs/Toronto_fdensity.png"))

## Fig. 10(b) spatial distribution signal relative l2 approximation error by various methods
DVEC = signal_transform_coeff(f, ht_elist_PC, ht_elist_VM, wavelet_packet_PC, wavelet_packet_VM, 𝚽, W, X)
approx_error_plot2(DVEC); approx_error_plt = plot!(legend = :topright, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
savefig(approx_error_plt, joinpath(@__DIR__, "../paperfigs/Toronto_fdensity_DAG_approx.png"))

## Fig. 11 spatial distribution signal 9 most important VM-NGWP vectors (ignore the DC vector)
parent_VM = HTree_findParent(ht_elist_VM); Wav_VM = best_basis_selection(f, wavelet_packet_VM, parent_VM); dvec_VM = Wav_VM' * f
importance_idx = sortperm(abs.(dvec_VM), rev = true)
for i = 2:10
    gplot(W, X; width=1); scatter_gplot!(X; marker = Wav_VM[:,importance_idx[i]], smallValFirst = true, ms = 3); important_NGW_basis_vectors = plot!(grid = false)
    savefig(important_NGW_basis_vectors, joinpath(@__DIR__, "../paperfigs/Toronto_fdensity_DAG_VM_NGW_important_basis_vector$(i).png"))
end

## Show some important HGLET vectors
G_Sig = GraphSig(1.0*W, xy=X, f=reshape(f, (N,1))); G_Sig = Adj2InvEuc(G_Sig);
GP = partition_tree_fiedler(G_Sig,:Lrw)
dmatrixH, dmatrixHrw, dmatrixHsym = HGLET_Analysis_All(G_Sig, GP) # expansion coefficients of 3-way HGLET bases
dvec_hglet, BS_hglet, trans_hglet = HGLET_GHWT_BestBasis(GP, dmatrixH = dmatrixH, dmatrixHrw = dmatrixHrw, dmatrixHsym = dmatrixHsym, costfun = 1) # best-basis among all combinations of bases
importance_idx = sortperm(dvec_hglet.^2; rev = true)[1:10]
for i = 2:10
    important_HGLET_basis_vector, _ = HGLET_Synthesis(reshape(spike(importance_idx[i],N), N, 1), GP, BS_hglet, G_Sig)
    w = important_HGLET_basis_vector[:,1]
    sgn = (maximum(w) > -minimum(w)) * 2 - 1
    gplot(W, X; width=1); scatter_gplot!(X; marker = sgn .* w, smallValFirst = true, ms = 3); important_HGLET_basis_vectors = plot!(grid = false)
    savefig(important_HGLET_basis_vectors, joinpath(@__DIR__, "../paperfigs/Toronto_fdensity_HGLET_important_basis_vector$(i).png"))
end

for i = 2:10
    println("(j, k, l) = ", HGLET_jkl(GP, BS_hglet.levlist[importance_idx[i]][1], BS_hglet.levlist[importance_idx[i]][2]))
end

## Fig. 12(a) pedestrian volume graph signal
fp = load(joinpath(@__DIR__, "..", "datasets", "new_toronto.jld"),"fp")
gplot(W, X; width=1); scatter_gplot!(X; marker = fp, smallValFirst = true, ms = 3); signal_plt = plot!(aspect_ratio = 1, grid = false)
savefig(signal_plt, joinpath(@__DIR__, "../paperfigs/Toronto_fp.png"))

## Fig. 12(b) pedestrian signal relative l2 approximation error by various methods
DVEC = signal_transform_coeff(fp, ht_elist_PC, ht_elist_VM, wavelet_packet_PC, wavelet_packet_VM, 𝚽, W, X)
approx_error_plot2(DVEC); approx_error_plt = plot!(legend = :topright, xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
savefig(approx_error_plt, joinpath(@__DIR__, "../paperfigs/Toronto_fp_DAG_approx.png"))

## Fig. 13 pedestrian signal 9 most important VM-NGWP vectors (ignore the DC vector)
parent_VM = HTree_findParent(ht_elist_VM); Wav_VM = best_basis_selection(fp, wavelet_packet_VM, parent_VM); dvec_VM = Wav_VM' * fp
importance_idx = sortperm(abs.(dvec_VM), rev = true)
for i = 2:10
    gplot(W, X; width=1); scatter_gplot!(X; marker = Wav_VM[:,importance_idx[i]], smallValFirst = true, ms = 3); important_NGW_basis_vectors = plot!(grid = false)
    savefig(important_NGW_basis_vectors, joinpath(@__DIR__, "../paperfigs/Toronto_fp_DAG_VM_NGW_important_basis_vector$(i).png"))
end

parent_PC = HTree_findParent(ht_elist_PC); Wav_PC = best_basis_selection(fp, wavelet_packet_PC, parent_PC); dvec_PC = Wav_PC' * fp
importance_idx = sortperm(abs.(dvec_PC), rev = true)
for i = 2:10
    w = Wav_PC[:,importance_idx[i]]
    sgn = (maximum(w) > -minimum(w)) * 2 - 1
    gplot(W, X; width=1); scatter_gplot!(X; marker = sgn .* w, smallValFirst = true, ms = 3); important_NGW_basis_vectors = plot!(grid = false)
    savefig(important_NGW_basis_vectors, joinpath(@__DIR__, "../paperfigs/Toronto_fp_DAG_PC_NGW_important_basis_vector$(i).png"))
end
