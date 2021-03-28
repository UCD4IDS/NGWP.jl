function shglet(𝚽, W, GP; ϵ = 0.3)
    rs = GP.rs
    inds = GP.inds
    (N, jmax) = Base.size(inds)

    GP.tag = zeros(Int, N, jmax)
    GP.tag[:, 1] = Vector{Int}(0:(N - 1))

    Uf = Matrix{Float64}(I, N, N)
    used_node = Set()

    HGLET = zeros(N, jmax, N)
    HGLET[:, 1, :] = 𝚽'

    sHGLET = zeros(N, jmax, N)
    sHGLET[:, 1, :] = 𝚽'
    for j = 2:jmax
        regioncount = count(!iszero, rs[:, j]) - 1
        # Uf = unitary_folding_operator(W, GP; ϵ = ϵ, J = j - 1)
        keep_folding!(Uf, used_node, W, GP; ϵ = ϵ, j = j - 1)
        for r = 1:regioncount
            indr = rs[r, j]:(rs[r + 1, j] - 1)
            GP.tag[indr, j] = Vector{Int}(0:(length(indr) - 1))
            HGLET[indr, j, :] = const_hglet_jk(W; idx = inds[indr, j])'
            sHGLET[indr, j, :] = const_shglet_jk(HGLET[indr, j, :]', Uf; idx = inds[indr, j])'
        end
    end
    return HGLET, sHGLET
end

function const_hglet_jk(W; idx = 1:size(Uf, 1))
    N = size(W, 1)
    Ljk = diagm(sum(W[idx, idx], dims = 1)[:]) - W[idx, idx]
    𝚽barjk = eigen(Matrix(Ljk)).vectors
    standardize_eigenvectors!(𝚽barjk)
    𝚽jk = zeros(N, length(idx))
    𝚽jk[idx, :] = 𝚽barjk
    return 𝚽jk
end



function const_shglet_jk(H, Uf; idx = 1:size(Uf, 1))
    # assemble smooth orthogonal projector
    P = Uf[idx, :]' * Uf[idx, :]
    return P * H
end
