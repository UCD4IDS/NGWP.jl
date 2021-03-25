println("Loading NGWP...")
using NGWP
using Test, JLD, MAT, Plots, LightGraphs, LinearAlgebra, SparseArrays, MTSG

#####################################################
# 1. Testing PC-NGWP and VM-NGWP on sunflower graph #
#####################################################
println("Testing NGWP on sunflower barbara eye signal")
## Build weighted sunflower graph for test
G, L, X = SunFlowerGraph(N = 400); N = nv(G)
𝛌, 𝚽 = eigen(Matrix(L)); standardize_eigenvectors!(𝚽)
Q = incidence_matrix(G; oriented = true)
W = 1.0 * adjacency_matrix(G)
edge_weight = [e.weight for e in edges(G)]

## Build Dual Graph by DAG metric
distDAG = eigDAG_Distance(𝚽, Q, N; edge_weight = edge_weight)
Gstar_Sig = dualgraph(distDAG)
G_Sig = GraphSig(W, xy = X)
GP_dual = partition_tree_fiedler(Gstar_Sig; swapRegion = false)
GP_primal = pairclustering(𝚽, GP_dual);

VM_NGWP = vm_ngwp(𝚽, GP_dual)
PC_NGWP = pc_ngwp(𝚽, GP_dual, GP_primal)

## use barbara eye signal for testing
using MAT
f = matread(joinpath(@__DIR__, "datasets",
                "sunflower_barbara_voronoi.mat"))["f_eye_voronoi"]
G_Sig.f = reshape(f, (N, 1))

DVEC = getall_expansioncoeffs(G_Sig, GP_dual, VM_NGWP, PC_NGWP, 𝚽)
DVEC_truth = JLD.load(joinpath(@__DIR__,
                        "datasets", "sunflower_barbara_feye_DVEC.jld"), "DVEC")

@testset "PC-NGWP Test" begin
    @test norm(DVEC[end-1] - DVEC_truth[end - 1]) / norm(DVEC_truth[end - 1]) < 1e-8
end

@testset "VM-NGWP Test" begin
    @test norm(DVEC[end] - DVEC_truth[end]) / norm(DVEC_truth[end]) < 1e-8
end
