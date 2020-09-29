"""
    eigTSD_Distance(Ve,V,lambda,Q,L;m = "Inf",dt = 0.1,tol = 1e-5)

EIGTSD\\_DISTANCE computes the TSDM distance matrix of Ve's column vectors on a graph.

# Input Argument
- `Ve::Matrix{Float64}`: feature matrix of eigenvectors, e.g., V.^2 or exp.(V)./sum(exp.(V), dims=1).
- `V::Matrix{Float64}`: matrix of graph Laplacian eigenvectors.
- `lambda::Array{Float64}`: vector of eigenvalues.
- `Q::Matrix{Float64}`: the oriented incidence matrix of the graph.
- `L::Matrix{Float64}`: the graph Laplacian matrix.
- `m*dt::Float64`: default is T = ∞, the stopping time T = m⋅dt in TSDM.
- `tol::Float64`: tolerance for convergence.

# Output Argument
- `dis::Matrix{Float64}`: distance matrix, d(φᵢ,φⱼ;T).

"""
function eigTSD_Distance(Ve,V,lambda,Q,L;m = "Inf",dt = 0.1,tol = 1e-5)
    N, J = size(Ve)
    dis = zeros(J,J)
    if m == "Inf"
        for i = 1:J-1, j = i+1:J
            cost = 0
            f₀ = Ve[:,i] - Ve[:,j]
            f = f₀
            c = V'*f₀
            global ind
            ind = 0
            while(norm(L*f,1)>tol)
                ind += 1
                cost += dt * norm(Q' * f,1)
                f = u_sol(c,V,lambda,ind*dt)
            end
            dis[i,j] = cost
        end
    else
        for i = 1:J-1, j = i+1:J
            cost = 0
            f₀ = Ve[:,i] - Ve[:,j]
            f = f₀
            c = V'*f₀
            for k = 1:m
                cost = cost + dt * norm(Q' * f,1)
                f = u_sol(c,V,lambda,k*dt)
            end
            dis[i,j] = cost
        end
    end
    return dis + dis'
end

"""
    TSD_Distance(p, q, V, lambda, Q, L; m = "Inf", dt = 0.1, tol = 1e-5)

TSD\\_DISTANCE computes the TSD distance between two vector meassures p and q on a graph.

# Input Argument
- `p::Array{Float64}`: the source vector measure.
- `q::Array{Float64}`: the destination vector measure.
- `V::Matrix{Float64}`: matrix of graph Laplacian eigenvectors.
- `lambda::Array{Float64}`: vector of eigenvalues.
- `Q::Matrix{Float64}`: the oriented incidence matrix of the graph.
- `L::Matrix{Float64}`: the graph Laplacian matrix.
- `m*dt::Float64`: default is T = ∞, the stopping time T = m⋅dt in TSDM.
- `tol::Float64`: tolerance for convergence.

# Output Argument
- `cost::Float64`: TSD distance d\\_TSD(p,q;T).

"""
function TSD_Distance(p,q,V,lambda,Q,L;m = "Inf",dt = 0.1,tol = 1e-5)
    cost = 0
    u₀ = q - p
    A = Q' * V
    b = V' * u₀
    if m == "Inf"
        t = 0
        while(true)
            t += dt
            increment = dt * norm(∇u(t, A, b, lambda), 1)
            cost += increment
            if increment < tol
                break
            end
        end
    else
        for k = 1:m
            cost += dt * norm(∇u(k*dt, A, b, lambda), 1)
        end
    end
    return cost
end

function ∇u(t, A, b, lambda)
    gu = A * (exp.(-t .* lambda) .* b)
    return gu
end

function u_sol(c,V,lambda,t)
    u = V * (exp.(-t .* lambda) .* c)
    return u
end
