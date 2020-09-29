using Optim

"""
    eigHAD_Distance(𝛷,lamb,numEigs)

EIGHAD\\_DISTANCE compute HAD "distance" (not really a distance) between pairwise graph Laplacian eigenvectors, i.e., d_HAD(𝜙ᵢ₋₁, 𝜙ⱼ₋₁) = 1/a_HAD(𝜙ᵢ₋₁, 𝜙ⱼ₋₁).

# Input Arguments
- `𝛷::Matrix{Float64}`: matrix of graph Laplacian eigenvectors, 𝜙ⱼ₋₁ (j = 1,...,size(𝛷,1)).
- `lamb::Array{Float64}`: array of eigenvalues. (ascending order)
- `indexEigs::Int`: default is all eigenvectors, indices of eigenvectors considered.

# Output Argument
- `dis::Matrix{Float64}`: a numEigs x numEigs affinity matrix, dis[i,j] = d_HAD(𝜙ᵢ₋₁, 𝜙ⱼ₋₁).
"""
function eigHAD_Distance(𝛷, lamb; indexEigs = 1:size(𝛷,2))
    A = eigHAD_Affinity(𝛷, lamb; indexEigs = indexEigs)
    n = size(A,1)
    dis = zeros(n,n)
    for i = 1:n, j = 1:n
        if A[i,j] == 0
            dis[i,j] = 1e9
        else
            dis[i,j] = 1/A[i,j]
        end
    end
    return dis
end

function eigHAD_Distance_neglog(𝛷, lamb; indexEigs = 1:size(𝛷,2))
    A = eigHAD_Affinity(𝛷, lamb; indexEigs = indexEigs)
    n = size(A,1)
    dis = zeros(n,n)
    for i = 1:n, j = 1:n
        if A[i,j] == 0
            dis[i,j] = 1e9
        else
            dis[i,j] = -log(A[i,j])
        end
    end
    return dis
end

"""
    eigHAD_Affinity(𝛷,lamb,numEigs)

EIGHAD_AFFINITY compute Hadamard (HAD) affinity between pairwise graph Laplacian eigenvectors.

# Input Arguments
- `𝛷::Matrix{Float64}`: matrix of graph Laplacian eigenvectors, 𝜙ⱼ₋₁ (j = 1,...,size(𝛷,1)).
- `lamb::Array{Float64}`: array of eigenvalues. (ascending order)
- `indexEigs::Int`: default is all eigenvectors, indices of eigenvectors considered.

# Output Argument
- `A::Matrix{Float64}`: a numEigs x numEigs affinity matrix, A[i,j] = a_HAD(𝜙ᵢ₋₁, 𝜙ⱼ₋₁).
"""
function eigHAD_Affinity(𝛷, lamb; indexEigs = 1:size(𝛷,2))
    N, numEigs = size(𝛷,1), length(indexEigs)
    indNoDC = setdiff(indexEigs, 1) # get rid of DC component
    J = length(indNoDC)
    A = zeros(J, J)
    for a in 1:J, b in a:J
        i, j = indNoDC[a], indNoDC[b]
        hadamardProd = 𝛷[:,i] .* 𝛷[:,j]
        if norm(hadamardProd,2) < 0.01/sqrt(N)
            continue
        end
        λ, μ = lamb[i], lamb[j]
        x₀ = 1 ./ (max(λ, μ))
        # Find minimizer t
        result = optimize(t -> abs(exp(-t[1]*λ) + exp(-t[1]*μ) - 1), [x₀], BFGS());
        t = Optim.minimizer(result)[1]
        # Compute Hadamard affinity
        heatEvolution = 𝛷 * Diagonal(exp.(-t .* lamb)) * 𝛷' * hadamardProd
        A[a,b] = norm(heatEvolution,2) / (norm(hadamardProd,2) + 1e-6)
    end
    A = A + A'; for i in 1:J; A[i,i] /= 2; end

    if 1 in indexEigs
        # Set affinity measure of 𝜙₀ with itself to be the maximum and equals to 1.
        A = matrix_augmentation(A)
        A[1,1] = maximum(A)
    end
    return A ./ maximum(A)
end

function matrix_augmentation(A)
    m, n = size(A)
    B = zeros(m+1, n+1)
    B[2:end, 2:end] = A
    return B
end
