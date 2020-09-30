## Set up environment and load required packages
using Pkg
Pkg.activate(".")

using NGWP, JLD, MAT, Plots, LightGraphs, MTSG
gr(dpi = 200)

## Build weighted sunflower graph
G, L, X = SunFlowerGraph(N = 400); N = nv(G)
lamb, 𝛷 = eigen(Matrix(L)); sgn = (maximum(𝛷, dims = 1)[:] .> -minimum(𝛷, dims = 1)[:]) .* 2 .- 1; 𝛷 = Matrix((𝛷' .* sgn)')
Q = incidence_matrix(G; oriented = true)
W = 1.0 * adjacency_matrix(G)
edge_weight = [e.weight for e in edges(G)]

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝛷,Q,N; edge_weight = edge_weight)
W_dual = sparse(dualGraph(distDAG)) #required: sparse dual weighted adjacence matrix

## Assemble natural graph wavelet packets
ht_elist_dual, ht_vlist_dual = HTree_EVlist(𝛷,W_dual)
wavelet_packet_dual = HTree_wavelet_packet(𝛷,ht_vlist_dual,ht_elist_dual)
ht_elist_varimax = ht_elist_dual
wavelet_packet_varimax = HTree_wavelet_packet_varimax(𝛷,ht_elist_varimax)

## Fig. 4(b) barbara eye graph signal
f = matread(joinpath(@__DIR__, "..", "datasets", "sunflower_barbara.mat"))["f_eye"]
scatter_gplot(X; marker = f, ms = LinRange(4.0, 14.0, N), smallValFirst = false, c = :greys); signal_plt = plot!(xlim = [-1.2,1.2], ylim = [-1.2,1.2], yflip = true, frame = :none)
savefig(signal_plt, "test/paperfigs/SunFlower_barbara_feye.png")

## Fig. 4(c) barbara eye relative l2 approximation error by various methods
DVEC = signal_transform_coeff(f, ht_elist_dual, ht_elist_varimax, wavelet_packet_dual, wavelet_packet_varimax, 𝛷, W, X)
approx_error_plot2(DVEC); plot!(legend = :topright); approx_error_plt = current()
savefig(approx_error_plt, "test/paperfigs/SunFlower_barbara_feye_DAG_approx.png")

## Fig. 5 barbara eye 9 most important VM-NGWP vectors (ignore the DC vector)
parent_varimax = HTree_findParent(ht_elist_varimax); Wav_varimax = best_basis_selection(f, wavelet_packet_varimax, parent_varimax); dvec_varimax = Wav_varimax' * f
importance_idx = sortperm(abs.(dvec_varimax), rev = true)
for i = 2:10
    scatter_gplot(X; marker = Wav_varimax[:,importance_idx[i]], ms = LinRange(4.0, 14.0, N), smallValFirst = false, c = :greys); important_NGW_basis_vectors = plot!(xlim = [-1.2,1.2], ylim = [-1.2,1.2], yflip = true, frame = :none)
    savefig(important_NGW_basis_vectors, "test/paperfigs/SunFlower_barbara_feye_DAG_VM_NGW_important_basis_vector$(i).png")
end

## Fig. 6(b) barbara pant graph signal
f = matread(joinpath(@__DIR__, "..", "datasets", "sunflower_barbara.mat"))["f_trouser"]
scatter_gplot(X; marker = f, ms = LinRange(4.0, 14.0, N), smallValFirst = false, c = :greys); signal_plt = plot!(xlim = [-1.2,1.2], ylim = [-1.2,1.2], yflip = true, frame = :none)
savefig(signal_plt, "test/paperfigs/SunFlower_barbara_ftrouser.png")

## Fig. 6(c) barbara eye relative l2 approximation error by various methods
DVEC = signal_transform_coeff(f, ht_elist_dual, ht_elist_varimax, wavelet_packet_dual, wavelet_packet_varimax, 𝛷, W, X)
approx_error_plot2(DVEC); plot!(legend = :topright); approx_error_plt = current()
savefig(approx_error_plt, "test/paperfigs/SunFlower_barbara_ftrouser_DAG_approx.png")

## Fig. 7 barbara eye 9 most important PC-NGWP vectors (ignore the DC vector)
parent_dual = HTree_findParent(ht_elist_dual); Wav_dual = best_basis_selection(f, wavelet_packet_dual, parent_dual); dvec_spectral = Wav_dual' * f
importance_idx = sortperm(abs.(dvec_spectral), rev = true)
for i = 2:10
    w = Wav_dual[:,importance_idx[i]]
    lm_sgn = (maximum(w) > -minimum(w)) * 2 - 1
    scatter_gplot(X; marker = lm_sgn .* w, ms = LinRange(4.0, 14.0, N), smallValFirst = false, c = :greys); important_NGW_basis_vectors = plot!(xlim = [-1.2,1.2], ylim = [-1.2,1.2], yflip = true, frame = :none)
    savefig(important_NGW_basis_vectors, "test/paperfigs/SunFlower_barbara_ftrouser_DAG_PC_NGW_important_basis_vector$(i).png")
end
