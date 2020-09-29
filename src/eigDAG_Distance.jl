"""
    eigDAG_Distance(V, Q, numEigs; edge_weights = 1)

EIGDAG_DISTANCE compute DAG distances between pairwise graph Laplacian eigenvectors.

# Input Arguments
- `V::Matrix{Float64}`: matrix of graph Laplacian eigenvectors, 𝜙ⱼ₋₁ (j = 1,...,size(V,1)).
- `Q::Matrix{Float64}`: incidence matrix of the graph.
- `numEigs::Int64`: number of eigenvectors considered.
- `edge_weight::Array{Float64}`: default value is 1, stands for unweighted graph (i.e., all edge weights equal to 1). For weighted graph, edge_weight is the weights vector, which stores the affinity weight of each edge.

# Output Argument
- `dis::Matrix{Float64}`: a numEigs x numEigs distance matrix, dis[i,j] = d_DAG(𝜙ᵢ₋₁, 𝜙ⱼ₋₁).

"""
function eigDAG_Distance(V,Q,numEigs; edge_weight = 1)
    dis = zeros(numEigs,numEigs)
    abs_∇V = abs.(Q' * V)
    for i = 1:numEigs, j = i+1:numEigs
        dis[i,j] = (edge_weight == 1) ? norm(abs_∇V[:,i]-abs_∇V[:,j],2) : sqrt(sum((abs_∇V[:,i]-abs_∇V[:,j]).^2 .* sqrt.(edge_weight)))
    end
    return dis + dis'
end

function eigDAG_Distance_normalized(V,Q,numEigs; edge_weight = 1)
    dis = zeros(numEigs,numEigs)
    abs_∇V = abs.(Q' * V)
    for i = 1:numEigs, j = i+1:numEigs
        dis[i,j] = (edge_weight == 1) ? norm(abs_∇V[:,i]-abs_∇V[:,j],2) : sqrt(sum((abs_∇V[:,i]-abs_∇V[:,j]).^2 .* sqrt.(edge_weight)))
        dis[i,j] /= norm(V[:,i] .* V[:,j],2)
    end
    return dis + dis'
end
