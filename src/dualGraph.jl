"""
    dualGraph(D; method = "inverse", σ = 1.0)

DUALGRAPH build the dual graph's weight matrix based on the given non-trivial eigenvector metric.

# Input Arguments
- `D::Matrix{Float64}`: eigenvector distance matrix
- `method::String`: default is by taking inverse of the distance between nodes. ways to build the dual graph edge weights. Option: inverse, Gaussian
- `σ::Float64`: default is 1. Gaussian parameter.

# Output Argument
- `W::Matrix{Float64}`: weight matrix of the dual graph

"""
function dualGraph(D; method = "inverse", σ = 1.0)
    n = size(D)[1]
    W = zeros(n,n)
    if method == "Gaussian"
        for i = 1:n-1, j = i+1:n
            W[i,j] = exp(- D[i,j] / σ^2)
        end
    elseif method == "inverse"
        for i = 1:n-1, j = i+1:n
            W[i,j] = 1/D[i,j]
        end
    end
    return W + W'
end
