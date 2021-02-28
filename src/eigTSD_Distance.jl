"""
    eigTSD_Distance(P::Matrix{Float64}, 𝚽::Matrix{Float64}, 𝛌::Vector{Float64},
                        Q::Matrix{Float64}, L::Matrix{Float64}; length::Any = 1,
                        T::Any = :Inf, dt::Float64 = 0.1, tol::Float64 = 1e-5)

computes the TSD distance matrix of P's column vectors on a graph.

# Input Argument
- `P::Matrix{Float64}`: vector measures, e.g., `𝚽.^2`.
- `𝚽::Matrix{Float64}`: matrix of the unweighted graph Laplacian eigenvectors.
- `𝛌::Vector{Float64}`: vector of eigenvalues.
- `Q::SparseMatrixCSC{Int64,Int64}`: the unweighted incidence matrix.
- `L::Matrix{Float64}`: the unweighted graph Laplacian matrix.
- `length::Any`: vector of edge lengths (default: `1` represents unweighted graphs)
- `T::Any`: the stopping time T in K_functional (default: `:Inf`)
- `dt::Float64`: time increment (default: `0.1`)
- `tol::Float64`: tolerance for convergence (default: `1e-5`)

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, d_{TSD}(φᵢ, φⱼ; T).

"""
function eigTSD_Distance(P::Matrix{Float64}, 𝚽::Matrix{Float64}, 𝛌::Vector{Float64},
                            Q::SparseMatrixCSC{Int64,Int64}, L::Matrix{Int}; length::Any = 1,
                            T::Any = :Inf, dt::Float64 = 0.1, tol::Float64 = 1e-5)
    N, ncols = Base.size(P)
    total_mass = sum(P, dims = 1)[:]
    if norm(total_mass - total_mass[1] * ones(ncols), Inf) > 10^4 * eps()
        @error("P's column measures do not share the same total mass.")
        return
    end
    # initialize the distance matrix
    dis = zeros(ncols, ncols)
    # store gradient of 𝚽 to avoid repeated computation
    ∇𝚽 = Q' * 𝚽

    for i = 1:(ncols - 1), j = (i + 1):ncols
        dis[i, j] = K_functional(P[:, i], P[:, j], ∇𝚽, 𝛌, L; length = length,
                                    T = T, dt = dt, tol = tol)[1]
    end
    return dis + dis'
end

"""
    K_functional(𝐩, 𝐪, 𝚽, 𝛌, Q, L; m = :Inf, dt = 0.1, tol = 1e-5)

computes the K_functional between two vector meassures 𝐩 and 𝐪 on a graph.

# Input Argument
- `𝐩::Vector{Float64}`: the source vector measure.
- `𝐪::Vector{Float64}`: the destination vector measure.
- `∇𝚽::Matrix{Float64}`: gradient of unweighted graph Laplacian eigenvectors.
- `𝛌::Vector{Float64}`: vector of eigenvalues.
- `L::Matrix{Int}`: the unweighted graph Laplacian matrix.
- `length::Any`: vector of edge lengths (default: 1 represents unweighted graphs)
- `T::Any`: the stopping time T in K_functional (default: :Inf)
- `dt::Float64`: time increment (default: 0.1)
- `tol::Float64`: tolerance for convergence (default: 1e-5)

# Output Argument
- `K::Float64`: TSD distance d_{TSD}(p, q; T).
- `t::Float64`: the actual stopping time

"""
function K_functional(𝐩::Vector{Float64}, 𝐪::Vector{Float64}, ∇𝚽::Matrix{Float64},
                        𝛌::Vector{Float64}, L::Matrix{Int}; length::Any = 1,
                        T::Any = :Inf, dt::Float64 = 0.1, tol::Float64 = 1e-5)
    if abs(sum(𝐩 - 𝐪)) > 10^4 * eps()
        @error("𝐩 and 𝐪 do not have the same total mass.")
    end
    K = 0
    f₀ = 𝐪 - 𝐩
    # store b to avoid repeated computation
    b = 𝚽' * f₀

    t = 0
    if T == :Inf
        # initialize `increment` as a very large number
        increment = BigFloat(typemax(Int32))
        while increment > tol
            t += dt
            increment = dt * weighted_1norm(∇f(t, ∇𝚽, b, 𝛌), length)
            K += increment
        end
    elseif isa(T, Float64)
        while t < T
            t += dt
            K += dt * weighted_1norm(∇f(t, ∇𝚽, b, 𝛌), length)
        end
    else
        @error("m can only be :Int or positive integers.")
    end
    return K, t
end

function ∇f(t, ∇𝚽, b, 𝛌)
    gu = ∇𝚽 * (exp.(-t * 𝛌) .* b)
    return gu
end

function f_sol(c, 𝚽, 𝛌, t)
    f = 𝚽 * (exp.(-t * 𝛌) .* c)
    return f
end

function weighted_1norm(x, length)
    if length == 1
        return norm(x, 1)
    else
        return sum(abs.(x) .* length)
    end
end
