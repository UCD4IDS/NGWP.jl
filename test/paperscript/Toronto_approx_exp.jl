using NGWP, JLD, MAT, Plots, LightGraphs, Distances, MTSG
gr(dpi = 200)

## Build weighted toronto street network graph
G = loadgraph(joinpath(@__DIR__, "..", "datasets", "new_toronto_graph.lgz")); N = nv(G)
X = load(joinpath(@__DIR__, "..", "datasets", "new_toronto.jld"),"xy")
W = 1.0 .* adjacency_matrix(G)
dist_X = pairwise(Euclidean(),X; dims = 1)
Weight = W .* dualGraph(dist_X; method = "inverse") # weighted adjacence matrix
L = Matrix(Diagonal(sum(Weight;dims = 1)[:]) - Weight)
lamb, 𝛷 = eigen(L); sgn = (maximum(𝛷, dims = 1)[:] .> -minimum(𝛷, dims = 1)[:]) .* 2 .- 1; 𝛷 = Matrix((𝛷' .* sgn)')
Q = incidence_matrix(G; oriented = true)
edge_weight = 1 ./ sqrt.(sum((Q' * X).^2, dims = 2)[:])

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝛷,Q,N; edge_weight = edge_weight)
W_dual = sparse(dualGraph(distDAG)) #required: sparse dual weighted adjacence matrix

## Assemble natural graph wavelet packets
ht_elist_dual, ht_vlist_dual = HTree_EVlist(𝛷,W_dual)
wavelet_packet_dual = HTree_wavelet_packet(𝛷,ht_vlist_dual,ht_elist_dual)
ht_elist_varimax = ht_elist_dual
wavelet_packet_varimax = HTree_wavelet_packet_varimax(𝛷,ht_elist_varimax)

## Fig. 8(a) a smooth spatial distribution of the street intersections graph signal
f = zeros(N); for i in 1:N; f[i] = length(findall(dist_X[:,i] .< 1/minimum(edge_weight))); end #fneighbor
gplot(W, X; width=1); scatter_gplot!(X; marker = f, smallValFirst = true, ms = 3); signal_plt = plot!(aspect_ratio = 1, grid = false)
savefig(signal_plt, "test/paperfigs/Toronto_fdensity.png")

## Fig. 8(b) spatial distribution signal relative l2 approximation error by various methods
DVEC = signal_transform_coeff(f, ht_elist_dual, ht_elist_varimax, wavelet_packet_dual, wavelet_packet_varimax, 𝛷, W, X)
approx_error_plot2(DVEC); plot!(legend = :topright); approx_error_plt = current()
savefig(approx_error_plt, "test/paperfigs/Toronto_fdensity_DAG_approx.png")

## Fig. 9 spatial distribution signal 9 most important VM-NGWP vectors (ignore the DC vector)
parent_varimax = HTree_findParent(ht_elist_varimax); Wav_varimax = best_basis_selection(f, wavelet_packet_varimax, parent_varimax); dvec_varimax = Wav_varimax' * f
importance_idx = sortperm(abs.(dvec_varimax), rev = true)
for i = 2:10
    gplot(W, X; width=1); scatter_gplot!(X; marker = Wav_varimax[:,importance_idx[i]], smallValFirst = true, ms = 3); important_NGW_basis_vectors = plot!(grid = false)
    savefig(important_NGW_basis_vectors, "test/paperfigs/Toronto_fdensity_DAG_VM_NGW_important_basis_vector$(i).png")
end

## Fig. 10(a) pedestrian volume graph signal
fp = load(joinpath(@__DIR__, "..", "datasets", "new_toronto.jld"),"fp")
gplot(W, X; width=1); scatter_gplot!(X; marker = fp, smallValFirst = true, ms = 3); signal_plt = plot!(aspect_ratio = 1, grid = false)
savefig(signal_plt, "test/paperfigs/Toronto_fp.png")

## Fig. 10(b) pedestrian signal relative l2 approximation error by various methods
DVEC = signal_transform_coeff(fp, ht_elist_dual, ht_elist_varimax, wavelet_packet_dual, wavelet_packet_varimax, 𝛷, W, X)
approx_error_plot2(DVEC); plot!(legend = :topright); approx_error_plt = current()
savefig(approx_error_plt, "test/paperfigs/Toronto_fp_DAG_approx.png")

## Fig. 11 barbara eye 9 most important PC-NGWP vectors (ignore the DC vector)
parent_varimax = HTree_findParent(ht_elist_varimax); Wav_varimax = best_basis_selection(fp, wavelet_packet_varimax, parent_varimax); dvec_varimax = Wav_varimax' * fp
importance_idx = sortperm(abs.(dvec_varimax), rev = true)
for i = 2:10
    gplot(W, X; width=1); scatter_gplot!(X; marker = Wav_varimax[:,importance_idx[i]], smallValFirst = true, ms = 3); important_NGW_basis_vectors = plot!(grid = false)
    savefig(important_NGW_basis_vectors, "test/paperfigs/Toronto_fp_DAG_VM_NGW_important_basis_vector$(i).png")
end

## Fig. 12 barbara eye 9 most important PC-NGWP vectors (ignore the DC vector)
parent_dual = HTree_findParent(ht_elist_dual); Wav_dual = best_basis_selection(fp, wavelet_packet_dual, parent_dual); dvec_pc = Wav_dual' * fp
importance_idx = sortperm(abs.(dvec_pc), rev = true)
for i = 2:10
    w = Wav_dual[:,importance_idx[i]]
    lm_sgn = (maximum(w) > -minimum(w)) * 2 - 1
    gplot(W, X; width=1); scatter_gplot!(X; marker = lm_sgn .* w, smallValFirst = true, ms = 3); important_NGW_basis_vectors = plot!(grid = false)
    savefig(important_NGW_basis_vectors, "test/paperfigs/Toronto_fp_DAG_PC_NGW_important_basis_vector$(i).png")
end
