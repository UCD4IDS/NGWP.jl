using NGWP
using Test, JLD, MAT, Plots, LightGraphs, LinearAlgebra, SparseArrays

#####################################################
# 1. Testing PC-NGWP and VM-NGWP on sunflower graph #
#####################################################
println("Testing NGWP on sunflower barbara eye signal")
## Build weighted sunflower graph for test
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

## use barbara eye signal for testing
f = matread(joinpath(@__DIR__, "datasets", "sunflower_barbara.mat"))["f_eye"]

DVEC = signal_transform_coeff(f, ht_elist_dual, ht_elist_varimax, wavelet_packet_dual, wavelet_packet_varimax, 𝛷, W, X)

DVEC_truth = JLD.load(joinpath(@__DIR__, "datasets", "sunflower_barbara_feye_DVEC.jld"), "DVEC")
dvec_PC_NGWP_truth = DVEC_truth[7]
dvec_VM_NGWP_truth = DVEC_truth[8]

@testset "PC-NGWP Test" begin
    @test norm(DVEC[7]-dvec_PC_NGWP_truth)/norm(dvec_PC_NGWP_truth) < 1e-8
end

@testset "VM-NGWP Test" begin
    @test norm(DVEC[8]-dvec_VM_NGWP_truth)/norm(dvec_VM_NGWP_truth) < 1e-8
end
