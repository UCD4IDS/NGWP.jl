# require: JuMP.jl and Clp.jl
using JuMP, Clp
"""
    eigROT_Distance(P, Q; edge_length = 1, α = 1.0)

EIGROT\\_DISTANCE computes the ROT distance matrix of P's column vectors on a graph.

# Input Argument
- `P::Matrix{Float64}`: a matrix whose columns are probability measures.
- `Q::Matrix{Float64}`: the oriented incidence matrix of the graph.
- `edge_length::Array{Float64}`: default is 1, which is for unweighted input graph. For weighted graph, edge_length is the length vector, which stores the length of each edge.
- `α::Float64`: default is 1.0. ROT parameter.

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, dis[i,j] = d\\_ROT(pᵢ, pⱼ; α).

"""
function eigROT_Distance(P, Q; edge_length = 1, α = 1.0)
    n = size(P,2)
    dis = zeros(n,n)
    le2 = [edge_length;edge_length]
    Q2 = [Q -Q]
    m2 = size(Q2,2)
    for i = 1:n-1, j = i+1:n
        f = P[:,i] - P[:,j]
        md = Model(optimizer_with_attributes(Clp.Optimizer, "LogLevel" => 0))
        @variable(md, w[1:m2] >= 0.0)
        edge_length == 1 ? @objective(md, Min, sum(w)) : @objective(md, Min, sum(w .* le2))
        @constraint(md, Q2 * w .== f)
        JuMP.optimize!(md)
        wt = abs.(JuMP.value.(w))
        dis[i,j] = edge_length == 1 ? norm(wt .^ α, 1) : norm((wt .^ α) .* le2, 1)
    end
    return dis + dis'
end


"""
    ROT_Distance(A, B, Q; edge_length = 1, α = 1.0)

ROT\\_DISTANCE computes the ROT distance matrix from A's column vectors to B's column vectors. If A, B are vector inputs, then it also returns the optimal transport plan solution and cost value.

# Input Argument
- `A::Matrix{Float64}`: a matrix whose columns are initial probability measures.
- `B::Matrix{Float64}`: a matrix whose columns are terminal probability measures.
- `Q::Matrix{Float64}`: the oriented incidence matrix of the graph.
- `edge_length::Array{Float64}`: default is 1, which is for unweighted input graph. For weighted graph, edge_length is the length vector, which stores the length of each edge.
- `α::Float64`: default is 1.0. ROT parameter.

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, dis[i,j] = d\\_ROT(aᵢ, bⱼ; α).

"""
function ROT_Distance(A, B, Q; edge_length = 1, α = 1.0)
    m = ndims(A) > 1 ? size(A,2) : 1
    n = ndims(B) > 1 ? size(B,2) : 1
    dis = zeros(m,n)
    le2 = [edge_length;edge_length]
    Q2 = [Q -Q]
    m2 = size(Q2,2)
    for i = 1:m, j = 1:n
        f = (ndims(A) > 1 && ndims(B) > 1) ? B[:,j] - A[:,i] : B - A
        md = Model(optimizer_with_attributes(Clp.Optimizer, "LogLevel" => 0))
        @variable(md, w[1:m2] >= 0.0)
        edge_length == 1 ? @objective(md, Min, sum(w)) : @objective(md, Min, sum(w .* le2))
        @constraint(md, Q2 * w .== f)
        JuMP.optimize!(md)
        wt = abs.(JuMP.value.(w))
        dis[i,j] = edge_length == 1 ? norm(wt .^ α, 1) : norm((wt .^ α) .* le2, 1)
        if ndims(A) == 1 && ndims(B) == 1
            return wt, dis[1,1]
        end
    end
    return dis
end


using OptimalTransport
"""
    eigEMD_Distance(P, C)

EIGEMD\\_DISTANCE computes the EMD distance matrix of P's column vectors on a graph.

# Input Argument
- `P::Matrix{Float64}`: a matrix whose columns are probability measures.
- `C::Matrix{Float64}`: the cost matrix, which is the ground metric between pairwise nodes.

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, dis[i,j] = d\\_EMD(pᵢ, pⱼ).

"""
function eigEMD_Distance(P, C)
    n = size(P,2)
    dis = zeros(n,n)
    for i = 1:n-1, j = i+1:n
        dis[i,j] = emd2(P[:,i], P[:,j], C)
    end
    return dis + dis'
end
