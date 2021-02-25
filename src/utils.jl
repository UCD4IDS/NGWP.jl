# project x onto matrix A's column space, given A is full column rank
function Proj(x,A)
    y = try
        pinv(A' * A) * (A' * x)
    catch err
        if err != Nothing
            A' * x
        end
    end
    return A * y
end

function wavelet_perp_Matrix(w,A)
    N, k = size(A)
    B = zeros(N, k - 1)
    if k == 1
        return
    end
    tmp = A'*w
    if tmp[1] > 1e-6
        M = Matrix{Float64}(I, k-1, k-1)
        return A*vcat((-tmp[2:end] ./ tmp[1])', M)
    end
    if tmp[end] > 1e-6
        M = Matrix{Float64}(I, k-1, k-1)
        return A*vcat(M,(-tmp[1:end-1] ./ tmp[end])')
    end
    ind = findall(tmp .== maximum(tmp))[1]
    rest_ind = setdiff([i for i in 1:k], ind)
    M = Matrix{Float64}(I, k-1, k-1)
    B = A * vcat(vcat(M[1:ind-1,:], (-tmp[rest_ind] ./ tmp[ind])'), M[ind:end,:])
    return B
end


# Gram-Schmidt Process Orthogonalization
function gram_schmidt(A; tol = 1e-12)
    # Input: matirx A
    # Output: orthogonalization matrix of A's column vectors

    # Convert the matrix to a list of column vectors
    a = [A[:,i] for i in 1:size(A,2)]

    # Start Gram-Schmidt process
    q = []
    complement_dim = 0
    for i = 1:length(a)
        qtilde = a[i]
        for j = 1:(i - 1)
            qtilde -= (q[j]' * a[i]) * q[j]
        end
        if norm(qtilde) < tol
            complement_dim = size(A,2) - i + 1
            break
        end
        push!(q, qtilde / norm(qtilde))
    end
    Q = zeros(size(A,1), length(q))
    for i = 1:length(q)
        Q[:,i] .= q[i]
    end
    return Q, complement_dim
end

#
"""
    mgslp(A::Matrix{Float64}; tol::Float64 = 1e-12, p::Float64 = 1.0)

Modified Gram-Schmidt Process Orthogonalization with ℓᵖ pivoting algorithm (MGSLp)

# Input Arguments
- `A::Matrix{Float64}`: whose column vectors are to be orthogonalized.

# Output Argument
- `A::Matrix{Float64}`: orthogonalization matrix of A's column vectors based on ℓᵖ pivoting.
"""
function mgslp(A::Matrix{Float64}; tol::Float64 = 1e-12, p::Float64 = 1.0)
    n = size(A, 2)
    # Convert the matrix to a list of column vectors
    a = [A[:, i] for i in 1:n]

    # Start modified Gram-Schmidt process
    q = []
    v = copy(a)
    vec_lp_norm = [norm(v[i], p) for i in 1:n]
    complement_dim = 0
    for i = 1:n
        # Pivoting based on minimum ℓᵖ-norm
        idx = findmin(vec_lp_norm)[2]
        v[i], v[idx + i - 1] = v[idx + i - 1], v[i]
        # Check the linear dependency
        if norm(v[i]) < tol
            complement_dim = n - i + 1
            break
        end
        qtilde = v[i] / norm(v[i], 2)
        vec_lp_norm = []
        for j = (i + 1):n
            v[j] .-= (qtilde' * v[j]) * qtilde
            push!(vec_lp_norm, norm(v[j], p))
        end
        push!(q, qtilde)
    end
    Q = zeros(size(A,1), length(q))
    for i = 1:length(q)
        Q[:,i] .= q[i]
    end
    return Q, complement_dim
end
