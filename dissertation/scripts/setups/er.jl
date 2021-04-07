Random.seed!(1234)
N, M = 50, 200
G = erdos_renyi(N, M)
# println("Is the ER graph connected: ", is_connected(G))

L = Matrix(laplacian_matrix(G))
𝛌, 𝚽 = eigen(L)
Q = incidence_matrix(G; oriented = true)
∇𝚽 = Q' * 𝚽
