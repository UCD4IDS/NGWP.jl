function eigsROT_Distance(P, W, X; α = 1.0)
    G = SimpleGraph(W)
    if ne(G) + 1 != nv(G) || !is_connected(G)
        @error("input graph is not a tree.")
    end
    d = degree(G)
    # find index of junction nodes
    ijc = findall(d .> 2)
    # cut the graph into several disconnected subgraphs
    Wc = deepcopy(W)
    for i in ijc; Wc[i, :] .= 0; Wc[:, i] .= 0; end
    # find the indices for each subgraph
    Ind = find_subgraph_inds(Wc)
    # low dimensional pmfs
    𝚯 = Ind' * P
    # the centroids of subgraphs
    Xs = Diagonal(1 ./ sum(Ind, dims = 1)[:]) * Ind' * X
    # build Gs, i.e., the graph of subgraphs
    Gs = Graph(size(Ind, 2))
    Ns = nv(Gs)
    ijcs = []
    for k = 1:Ns
        supportind = findall(Ind[:, k] .== 1)
        if issubset(supportind, ijc)
            push!(ijcs, k)
        end
    end
    for k in setdiff(1:Ns, ijcs)
        supportind = findall(Ind[:,k] .== 1)
        for i = 1:length(ijc)
            if sum(W[supportind, ijc[i]]) > 0
                add_edge!(Gs, Edge(k, ijcs[i]))
            end
        end
    end
    Ws = 1.0 * adjacency_matrix(Gs)
    # compute the sROT distance matrix
    Qs = incidence_matrix(Gs; oriented = true)
    dist_sROT = eigROT_Distance(𝚯, Qs; edge_length = 1, α = α)
    return dist_sROT, Ws, Xs, 𝚯
end


function find_subgraph_inds(Wc)
    N = size(Wc, 1)
    Dc = Diagonal(sum(Wc, dims=1)[:])
    Lc = Matrix(Dc - Wc)
    𝛌c, 𝚽c = eigen(Lc)
    p = findall(𝛌c .< 100 * eps())
    Uc = 𝚽c[:, p]
    a = round.(1e8 * Uc * Uc' * [k for k = 1:N])
    ind = unique(a)
    Ind = zeros(N, length(ind))
    for k = 1:length(ind); Ind[findall(a .== ind[k]), k] .= 1; end
    return Ind
end
