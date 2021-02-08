println("Loading NGWP...")
using NGWP
using Test, JLD, MAT, Plots, LightGraphs, LinearAlgebra, SparseArrays

#####################################################
# 1. Testing PC-NGWP and VM-NGWP on sunflower graph #
#####################################################
println("Testing NGWP on sunflower barbara eye signal")
## Build weighted sunflower graph for test
G, L, X = SunFlowerGraph(); N = nv(G)
𝛌, 𝚽 = eigen(Matrix(L)); sgn = (maximum(𝚽, dims = 1)[:] .> -minimum(𝚽, dims = 1)[:]) .* 2 .- 1; 𝚽 = 𝚽 * Diagonal(sgn);
Q = incidence_matrix(G; oriented = true)
W = 1.0 * adjacency_matrix(G)
edge_weight = [e.weight for e in edges(G)]

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝚽,Q,N; edge_weight = edge_weight)
W_dual = sparse(dualGraph(distDAG)) #required: sparse dual weighted adjacence matrix

## Assemble natural graph wavelet packets
ht_elist_PC, ht_vlist_PC = HTree_EVlist(𝚽,W_dual)
wavelet_packet_PC = HTree_wavelet_packet(𝚽,ht_vlist_PC,ht_elist_PC)
ht_elist_VM = ht_elist_PC
wavelet_packet_VM = HTree_wavelet_packet_varimax(𝚽,ht_elist_VM)

## use barbara eye signal for testing
f = matread(joinpath(@__DIR__, "datasets", "sunflower_barbara_voronoi.mat"))["f_eye_voronoi"]

DVEC = signal_transform_coeff(f, ht_elist_PC, ht_elist_VM, wavelet_packet_PC, wavelet_packet_VM, 𝚽, W, X)
DVEC_truth = JLD.load(joinpath(@__DIR__, "datasets", "sunflower_barbara_feye_DVEC.jld"), "DVEC")

@testset "PC-NGWP Test" begin
    @test norm(DVEC[end-1]-DVEC_truth[end-1])/norm(DVEC_truth[end-1]) < 1e-8
end

@testset "VM-NGWP Test" begin
    @test norm(DVEC[end]-DVEC_truth[end])/norm(DVEC_truth[end]) < 1e-8
end
