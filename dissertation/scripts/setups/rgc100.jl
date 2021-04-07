using Plots, LightGraphs, JLD, LaTeXStrings, MTSG, NGWP
using Plots.PlotMeasures

G = loadgraph("../datasets/RGC100.lgz"); N = nv(G)
X = load("../datasets/RGC100_xyz.jld", "xyz")[:, 1:2]
X3 = load("../datasets/RGC100_xyz.jld", "xyz")
L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L)
standardize_eigenvectors!(𝚽)
W = 1.0 * adjacency_matrix(G)
Q = incidence_matrix(G; oriented = true)
∇𝚽 = Q' * 𝚽
edge_length = sqrt.(sum((Q' * X3).^2, dims = 2)[:])
