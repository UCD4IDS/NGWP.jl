"""
    natural_eigdist(𝚽, 𝛌, Q; α = 1.0, T = :Inf, dt = 0.01,
                    input_format = :eigenvectors, distance = :DAG,
                    edge_weight = 1, edge_length = 1)

compute natural distances between graph Laplacian eigenvectors.

# Input Arguments
- `𝚽::Matrix{Float64}`: matrix of (weighted) graph Laplacian eigenvectors.
- `𝛌::Vector{Float64}`: vector of eigenvalues.
- `Q::Matrix{Float64}`: unweighted incidence matrix of the graph.
- `α::Float64`: ROT parameter. (default: `1.0`)
- `T::Any`: TSD parameter, i.e., the stopping time T in K_functional (default: `:Inf`)
- `dt::Float64`: TSD parameter, i.e., the time increment (default: `0.1`)
- `input_format::Symbol`: options: `:eigenvectors`, `:zero_measures`, `:pmf1`
    and `:pmf2` (default: `:eigenvectors`)
- `distance::Symbol`: options: `:ROT`, `:HAD`, `:DAG` and `:TSD` (default: `:DAG`)
- `edg_length::Any`: vector of edge lengths (default: 1 represents unweighted graphs)
- `edge_weight::Any`: the weights vector, which stores the affinity weight of
    each edge (default: `1` represents unweighted graphs).

# Output Argument
- `dis::Matrix{Float64}`: the distance matrix, dis[i,j] = d(𝜙ᵢ₋₁, 𝜙ⱼ₋₁).

"""
function natural_eigdist(𝚽, 𝛌, Q; α = 1.0, T = :Inf, dt = 0.01,
                         input_format = :eigenvectors, distance = :DAG,
                         edge_weight = 1, edge_length = 1)
    N = size(Q, 1)
    if input_format == :eigenvectors
        P = 𝚽
    elseif input_format == :zero_measures
        P = 𝚽
        P[:, 1] .= 0
    elseif input_format == :pmf1
        P = 𝚽.^2
    elseif input_format == :pmf2
        P = exp.(𝚽) ./ sum(exp.(𝚽), dims = 1)
    else
        @error("input_format does not support $(input_format)!")
        return
    end

    if distance == :ROT
        D = eigROT_Distance(P, Q; edge_length = edge_length, α = α)
    elseif distance == :HAD
        D = eigHAD_Distance(P, 𝛌)
    elseif distance == :DAG
        D = eigDAG_Distance(P, Q, N; edge_weight = edge_weight)
    elseif distance == :TSD
        L = Q * Q'
        t𝛌, t𝚽 = eigen(Matrix(L))
        D = eigTSD_Distance(P, t𝚽, t𝛌, Q; length = edge_length, T = T, dt = dt)
    else
        error("distance does not $(distance)!")
        return
    end

    return D
end
