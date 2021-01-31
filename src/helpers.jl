"""
    findminimum(v, n)

FINDMINIMUM finds the first n smallest elements' indices.

# Input Arguments
- `v::Array{Float64}`: the candidate values for selection.
- `n::Int`: number of smallest elements for consideration.

# Output Argument
- `idx::Array{Int}`: n smallest elements' indices.

"""
function findminimum(v, n)
    idx = sortperm(v)[1:n]
    return idx
end


"""
    spike(i,n)

SPIKE gives the n-dim spike vector with i-th element equals 1.

# Input Arguments
- `i::Int`: index for one.
- `n::Int`: dimension of the target spike vector.

# Output Argument
- `a::Array{Float64}`: the n-dim spike vector with i-th element equals 1.

"""
function spike(i,n)
    a = zeros(n)
    a[i] = 1
    return a
end

"""
    characteristic(list,n)

CHARACTERISTIC gives the characteristic function in n-dim vector space with values of index in list equal to 1.

# Input Arguments
- `list::Array{Int}`: list of indices.
- `n::Int`: dimension of the target vector.

# Output Argument
- `v::Array{Float64}`: the n-dim characteristic vector with values of index in list equal to 1.

"""
function characteristic(list,n)
    v = zeros(n)
    v[list] .= 1.0
    return v
end


"""
    heat_sol(f0,Φ,Σ,t)

HEAT\\_SOL gives the solution of heat partial differential equation with initial condition u(⋅, 0) = f0

# Input Arguments
- `f0::Array{Float64}`: initial condition vector.
- `Φ::Matrix{Float64}`: graph Laplacian eigenvectors, served as graph Fourier transform matrix
- `Σ::Array{Int}`: diagonal matrix of eigenvalues.
- `t::Float`: time elapse.

# Output Argument
- `u::Array{Float64}`: the solution vector at time t

"""
function heat_sol(f0,Φ,Σ,t)
    u = Φ * (exp.(-t .* Σ) .* Φ' * f0)
    return u
end


"""
    freq_band_matrix(ls,n)

FREQ\\_BAND\\_MATRIX provides characteristic diagonal matrix, which is useful for spectral graph filters design.

# Input Arguments
- `ls::Array{Int}`: list of indices.
- `n::Int`: dimension of the target vector.

# Output Argument
- `D::Array{Float64}`: the zero/one diagonal matrix.

"""
function freq_band_matrix(ls, n)
    f = characteristic(list,n)
    return Diagonal(f)
end


"""
    scatter_gplot(X; marker = nothing, ms = 4, plotOrder = :normal, c = :viridis)

SCATTER\\_GPLOT generates a scatter plot figure, which is for quick viewing of a graph signal.
SCATTER\\_GPLOT!(X; ...) adds a plot to `current` one.

# Input Arguments
- `X::Matrix{Float64}`: points locations, can be 2-dim or 3-dim.
- `marker::Array{Float64}`: default is nothing. Present different colors given different signal value at each node.
- `ms::Array{Float64}`: default is 4. Present different node sizes given different signal value at each node.
- `plotOrder::Symbol`: default is normal. Optional choices :s2l or :l2s, i.e., plots from the smallest value of `marker` to the largest value or the other way around.
- `c::Symbol`: default is :viridis. Colors.

"""
function scatter_gplot(X; marker = nothing, ms = 4, plotOrder = :normal, c = :viridis)
    dim = size(X,2)
    if marker != nothing && plotOrder != :normal
        if plotOrder == :s2l
            idx = sortperm(marker)
        elseif plotOrder == :l2s
            idx = sortperm(marker, rev=true)
        else
            print("Error: plotOrder only supports for :normal, :s2l, or :l2s.")
        end
        X = X[idx,:]
        marker = marker[idx]
        if length(ms) > 1
            ms = ms[idx]
        end
    end
    if dim == 2
        scatter(X[:,1],X[:,2], marker_z = marker, ms = ms, c = c, legend = false, mswidth = 0, cbar = true, aspect_ratio = 1, grid = false)
    elseif dim == 3
        scatter(X[:,1],X[:,2],X[:,3], marker_z = marker, ms = ms, c = c, legend = false, mswidth = 0, cbar = true, aspect_ratio = 1, grid = false)
    else
        print("Dimension Error: scatter_gplot only supports for 2-dim or 3-dim scatter plots.")
    end
end

function scatter_gplot!(X; marker = nothing, ms = 4, plotOrder = :normal, c = :viridis)
    dim = size(X,2)
    if marker != nothing && plotOrder != :normal
        if plotOrder == :s2l
            idx = sortperm(marker)
        elseif plotOrder == :l2s
            idx = sortperm(marker, rev=true)
        else
            print("Error: plotOrder only supports for :normal, :s2l, or :l2s.")
        end
        X = X[idx,:]
        marker = marker[idx]
        if length(ms) > 1
            ms = ms[idx]
        end
    end
    if dim == 2
        scatter!(X[:,1],X[:,2], marker_z = marker, ms = ms, c = c, legend = false, mswidth = 0, cbar = true, aspect_ratio = 1, grid = false)
    elseif dim == 3
        scatter!(X[:,1],X[:,2],X[:,3], marker_z = marker, ms = ms, c = c, legend = false, mswidth = 0, cbar = true, aspect_ratio = 1, grid = false)
    else
        print("Dimension Error: scatter_gplot! only supports for 2-dim or 3-dim scatter plots.")
    end
end



"""
    cat_plot(X; marker = nothing, ms = 4)

CAT\\_PLOT generates a scatter plot figure for cat example, which is for quick viewing of a graph signal within a specific range (i.e., xlims, ylims, zlims).
CAT\\_PLOT!(X; ...) adds a plot to `current` one.

# Input Arguments
- `X::Matrix{Float64}`: 3-dim points.
- `marker::Array{Float64}`: default is nothing. Present different colors given different signal value at each node.
- `ms::Array{Float64}`: default is 4. Present different node sizes given different signal value at each node.

"""
function cat_plot(X; marker = nothing, ms = 4)
    scatter(X[:,1],X[:,2],X[:,3], marker_z = marker, ms = ms, c = :viridis, legend = false, cbar = true, aspect_ratio = 1, xlims = [-100, 100], ylims = [-100, 100], zlims = [-100, 100])
end

function cat_plot!(X; marker = nothing, ms = 4)
    scatter!(X[:,1],X[:,2],X[:,3], marker_z = marker, ms = ms, c = :viridis, legend = false, cbar = true, aspect_ratio = 1, xlims = [-100, 100], ylims = [-100, 100], zlims = [-100, 100])
end


"""
    approx_error_plot(ortho_mx_list, f; fraction_cap = 0.3, label = false, Save = false, path = "")

APPROX\\_ERROR\\_PLOT draw approx. error figure w.r.t. fraction of kept coefficients

# Input Arguments
- `ortho_mx_list::Array{Matrix{Float64}}`: a list of orthonormal matrices.
- `f::Array{Float64}`: target graph signal for approximation.
- `fraction_cap::Float`: default is 0.3. The capital of fration of kept coefficients.

"""
function approx_error_plot(ortho_mx_list, f; fraction_cap = 0.3, label = false, Save = false, path = "")
    N = length(f)
    L = length(ortho_mx_list)
    err = [[1.0] for _ in 1:L]
    coeff = [mx'*f for mx in ortho_mx_list]

    for frac = 0.01:0.01:fraction_cap
        numKept = Int(ceil(frac * N))
        for l in 1:L
            ind = sortperm(coeff[l].^2, rev = true)[numKept+1:end]
            push!(err[l], norm(coeff[l][ind])/norm(f))
        end
    end

    gr(dpi = 300)
    fraction = 0:0.01:fraction_cap
    plt = plot(fraction, err, yscale=:log10, lab = label, linewidth = 3, xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    if Save
        savefig(plt, path)
        return "figure saved! @ " * path
    end
    return "use current() to show figure."
end

################################################################################
####################### Approximation error plot################################
################################################################################
### function to plot the approximation error curve
function approx_error_plot2(DVEC::Array{Array{Float64,1},1}; frac = 0.50)
    plot(xaxis = "Fraction of Coefficients Retained", yaxis = "Relative Approximation Error")
    T = ["Laplacian", "HGLET", "Haar", "Walsh",  "GHWT_c2f", "GHWT_f2c", "eGHWT", "PC-NGWP", "VM-NGWP"]
    L = [(:dashdot,:red), (:solid,:brown), (:dashdot,:orange), (:dashdot,:pink), (:solid,:gray), (:solid,:green), (:solid,:blue), (:solid,:purple), (:solid,:black)]
    for i = 1:length(DVEC)
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2) # the smallest first
        er = max.(sqrt.(reverse(cumsum(dvec_sort)))/dvec_norm, 1e-12) # this is the relative L^2 error of the whole thing, i.e., its length is N
        p = Int64(floor(frac*N)) + 1 # upper limit
        plot!(frac*(0:(p-1))/(p-1), er[1:p], yaxis=:log, xlims = (0.,frac), label = T[i], line = L[i], linewidth = 2, grid = false)
    end
end

function signal_transform_coeff(f, ht_elist_dual, ht_elist_varimax, wavelet_packet_dual, wavelet_packet_varimax, 𝛷, W, X)
    ## Best basis selection algorithm
    parent_dual = HTree_findParent(ht_elist_dual)
    Wav_dual = best_basis_selection(f, wavelet_packet_dual, parent_dual)
    parent_varimax = HTree_findParent(ht_elist_varimax)
    Wav_varimax = best_basis_selection(f, wavelet_packet_varimax, parent_varimax)
    ############# varimax NGW coefficients
    dvec_varimax = Wav_varimax' * f
    ############# spectral_prioritized PC NGW coefficients
    dvec_spectral = Wav_dual' * f
    ############# plain Laplacian eigenvectors coefficients
    dvec_Laplacian = 𝛷' * f
    ## MTSG
    G_Sig = GraphSig(1.0*W, xy=X, f=reshape(f, (N,1)))
    G_Sig = Adj2InvEuc(G_Sig)
    GP = partition_tree_fiedler(G_Sig,:Lrw)
    dmatrixH = HGLET_Analysis_All(G_Sig, GP)[1] # expansion coefficients of 3-way HGLET bases
    dvec_hglet, BS_hglet, trans_hglet = HGLET_GHWT_BestBasis(GP, dmatrixH = dmatrixH, costfun = 1) # best-basis among all combinations of bases

    dmatrix = ghwt_analysis!(G_Sig, GP=GP)
    ############# Haar
    BS_haar = bs_haar(GP)
    dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)
    ############# Walsh
    BS_walsh = bs_walsh(GP)
    dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)
    ############# GHWT_c2f
    dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)
    ############# GHWT_f2c
    dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)
    ############# eGHWT
    dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)
    DVEC = [dvec_Laplacian[:], dvec_hglet[:], dvec_haar[:], dvec_walsh[:], dvec_c2f[:], dvec_f2c[:], dvec_eghwt[:], dvec_spectral[:], dvec_varimax[:]]
    return DVEC
end

function signal_transform_coeff2(f, g, ht_elist_dual, ht_elist_varimax, wavelet_packet_dual, wavelet_packet_varimax, 𝛷, W, X)
    parent_dual = HTree_findParent(ht_elist_dual); W_PC = best_basis_selection(f, wavelet_packet_dual, parent_dual)
    parent_varimax = HTree_findParent(ht_elist_varimax); W_VM = best_basis_selection(f, wavelet_packet_varimax, parent_varimax)
    ############# spectral_prioritized PC NGW coefficients
    dvec_spectral = W_PC' * [f g]
    ############# varimax NGW coefficients
    dvec_varimax = W_VM' * [f g]
    ############# plain Laplacian eigenvectors coefficients
    dvec_Laplacian = 𝛷' * [f g]
    tmp=zeros(length(f),2); tmp[:,1]=f; tmp[:,2]=g; G_Sig=GraphSig(1.0*W, xy=X, f=tmp); G_Sig = Adj2InvEuc(G_Sig); GP = partition_tree_fiedler(G_Sig,:Lrw);
    dmatrix = ghwt_analysis!(G_Sig, GP=GP)
    ############# Haar
    BS_haar = bs_haar(GP)
    dvec_haar = dmatrix2dvec(dmatrix, GP, BS_haar)
    ############# Walsh
    BS_walsh = bs_walsh(GP)
    dvec_walsh = dmatrix2dvec(dmatrix, GP, BS_walsh)
    ############# GHWT_c2f
    dvec_c2f, BS_c2f = ghwt_c2f_bestbasis(dmatrix, GP)
    ############# GHWT_f2c
    dvec_f2c, BS_f2c = ghwt_f2c_bestbasis(dmatrix, GP)
    ############# eGHWT
    dvec_eghwt, BS_eghwt = ghwt_tf_bestbasis(dmatrix, GP)
    DVEC = [dvec_haar, dvec_walsh, dvec_Laplacian, dvec_c2f, dvec_f2c, dvec_eghwt, dvec_spectral, dvec_varimax]
    return DVEC
end

function integrate_approx_results(DVEC, num_kept_coeffs, filename)
    ERR = Array{Float64,1}[]
    for i in 1:length(DVEC)
        dvec = DVEC[i]
        N = length(dvec)
        dvec_norm = norm(dvec,2)
        dvec_sort = sort(dvec.^2)  # the smallest first
        er = max.(sqrt.(reverse(cumsum(dvec_sort))) ./ dvec_norm, 1e-12)  # this is the relative L^2-norm error of the whole thing, i.e., its length is N
        push!(ERR, er[num_kept_coeffs])
    end
    frames_approx_res = CSV.File(joinpath(@__DIR__, "..", "datasets", filename))
    er_soft_cluster_frame = [max(sqrt(frames_approx_res[i][2]), 1e-12) for i in 1:length(num_kept_coeffs)]
    push!(ERR, er_soft_cluster_frame)
    er_SGWT = [max(sqrt(frames_approx_res[i][3]), 1e-12) for i in 1:length(num_kept_coeffs)]
    push!(ERR, er_SGWT)
    return ERR
end

function integrate_approx_results2(DVEC, num_kept_coeffs, filename)
    ERR = Array{Float64,1}[]
    for i in 1:length(DVEC)
        dvec_f = DVEC[i][:,1]
        dvec_g = DVEC[i][:,2]
        N = length(dvec_f)
        dvec_f_norm = norm(dvec_f)
        ind = sortperm(dvec_g.^2, rev = true)  # the largest first
        er = zeros(N)
        for k in 1:N
            er[k] = max(sqrt(norm(dvec_f[ind[1:k]] - dvec_g[ind[1:k]])^2 + norm(dvec_f[ind[k+1:end]])^2)/dvec_f_norm, 1e-12)
        end
        push!(ERR, er[num_kept_coeffs])
    end
    frames_approx_res = CSV.File(joinpath(@__DIR__, "..", "datasets", filename))
    er_soft_cluster_frame = [max(sqrt(frames_approx_res[i][2]), 1e-12) for i in 1:length(frames_approx_res)]
    push!(ERR, er_soft_cluster_frame)
    er_SGWT = [max(sqrt(frames_approx_res[i][3]), 1e-12) for i in 1:length(frames_approx_res)]
    push!(ERR, er_SGWT)
    return ERR
end

function approx_error_plot3(ERR::Array{Array{Float64,1},1}; num_kept_coeffs = 10:10:280, LW = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2])
    gr(dpi = 400)
    plot(xaxis = "Number of Coefficients Retained", yaxis = "Relative Approximation Error")
    T = ["Haar", "Walsh", "Laplacian", "GHWT_c2f", "GHWT_f2c", "eGHWT", "PC-NGWP", "VM-NGWP", "SC-Frame", "SGWT"]
    L = [(:dashdot,:orange), (:dashdot,:pink), (:dashdot, :red), (:solid, :gray), (:solid, :green), (:solid, :blue), (:solid, :purple), (:solid, :black), (:dash, :navy), (:dash, :teal)]
    for i in 1:length(ERR)
        plot!(num_kept_coeffs, ERR[i], yaxis=:log, xlims = (0.,num_kept_coeffs[end]), label = T[i], line = L[i], linewidth = LW[i])
    end
end


"""
    sortWaveletsByCenteredLocations(Wav)

sort wavelets by centered locations

# Input Argument
- `Wav::Matrix{Float64}`: a matrix whose columns are wavelet vectors.

# Output Argument
- `Wav::Matrix{Float64}`: the sorted matrix.
"""
function sortWaveletsByCenteredLocations(Wav)
    ord = findmax(abs.(Wav), dims = 1)[2][:]
    idx = sortperm([i[1] for i in ord])
    return Wav[:,idx]
end

using Clustering
"""
    spectral_clustering(𝛷, M)

SPECTRAL_CLUSTERING return M graph clusters, i.e., {Vₖ| k = 1,2,...,M}.

# Input Argument
- `𝛷::Matrix{Float64}`: the matrix of graph Laplacian eigenvectors.
- `M::Int64`: the number of graph clusters.

# Output Argument
- `clusters::Array{Array{Int64}}`: graph cluster indices.

"""
function spectral_clustering(𝛷, M)
    if M < 2
        return [1:size(𝛷,1)]
    end
    cluster_indices = assignments(kmeans(𝛷[:,2:M]', M))
    clusters = Array{Array{Int64,1},1}()
    for k in 1:M
        push!(clusters, findall(cluster_indices .== k)[:])
    end
    return clusters
end

"""
    transform2D(X; s = 1, t = [0,0])

TRANSFORM2D dilate each point of `X` by scale s and translate by 2D vector t.
"""
function transform2D(X; s = 1, t = [0,0])
    X1 = X .* s
    X2 = zeros(size(X))
    for i in 1:size(X,1)
        X2[i,1] = X1[i,1] + t[1]
        X2[i,2] = X1[i,2] + t[2]
    end
    return X2
end

"""
    NN_rendering(X, Img_Mat)

NN\\_RENDERING generates a rendering signal at each point of `X` from the image `Img_Mat` by nearest neighbor method.
"""
function NN_rendering(X, Img_Mat)
    N = size(X,1)
    f = zeros(N)
    for i in 1:N
        nn_x, nn_y = Int(round(X[i, 2])), Int(round(X[i, 1]))
        if nn_x < 1 || nn_x > size(Img_Mat, 2) || nn_y < 1 || nn_y > size(Img_Mat, 1)
            print("Error: pixel out of boundary!")
            return
        end
        f[i] = Img_Mat[nn_x, nn_y]
    end
    return f
end

"""
    Bilinear_rendering(X, Img_Mat)

NN\\_RENDERING generates a rendering signal at each point of `X` from the image `Img_Mat` by bilinear interpolation method.
"""
function Bilinear_rendering(X, Img_Mat)
    N = size(X,1)
    f = zeros(N)
    for i in 1:N
        x1, x2, y1, y2 = Int(floor(X[i, 2])), Int(floor(X[i, 2])) + 1, Int(floor(X[i, 1])), Int(floor(X[i, 1])) + 1
        x, y = X[i,2], X[i,1]
        F = [Img_Mat[x1,y1] Img_Mat[x1,y2]
            Img_Mat[x2,y1] Img_Mat[x2,y2]]
        prod_res = 1/((x2 - x1) * (y2 - y1)) * [x2-x x-x1] * F * [y2-y y-y1]'
        f[i] = prod_res[1,1]
    end
    return f
end

"""
    dct1d(k, N)

DCT1D returns k-th 1D DCT basis vector in Rᴺ.

# Input Arguments
- `k::Int64`: ord of DCT basis vector. k = 1,2,...,N.
- `N::Int64`: vector dimension.

# Output Argument
- `φ::Array{Float64}`: k-th 1D DCT basis vector in Rᴺ. (k is 1-indexed)
"""
function dct1d(k, N)
    φ = [cos(π*(k-1)*(l+0.5)/N) for l = 0:N-1]
    return φ ./ norm(φ, 2)
end

"""
    dct2d_basis(N1, N2)

DCT2D\\_BASIS returns 2D DCT basis vectors in [0,1] x [0,1] with N1-1 and N2-1 subintervals respectively.

# Input Arguments
- `N1::Int64`: number of nodes in x-axis.
- `N2::Int64`: number of nodes in y-axis.

# Output Argument
- `𝚽::Matrix{Float64}`: 2D DCT basis vectors.
"""
function dct2d_basis(N1, N2)
    N = N1 * N2
    𝚽 = zeros(N, N)
    ind = 1
    for i in 1:N1, j in 1:N2
        φ₁, φ₂ = dct1d(i, N1), dct1d(j, N2)
        φ = reshape(φ₁*φ₂', N)
        𝚽[:,ind] = φ
        ind += 1
    end
    return 𝚽
end

"""
    alternating_numbers(n)

ALTERNATING\\_NUMBERS e.g., n = 5, returns [1,5,2,4,3]; n = 6, returns [1,6,2,5,3,4]

# Input Arguments
- `N1::Int64`: number of nodes in x-axis.

# Output Argument
- `arr::Array{Int64}`: result array.
"""
function alternating_numbers(n)
    mid = Int(ceil(n/2))
    arr1 = 1:mid
    arr2 = n:-1:(mid+1)
    arr = Array{Int64}(zeros(n))
    p1, p2 = 1, 1
    for i = 1:n
        if i % 2 == 1
            arr[i] = arr1[p1]
            p1 += 1
        else
            arr[i] = arr2[p2]
            p2 += 1
        end
    end
    return arr
end

"""
    compute_SNR(f, g)

COMPUTE\\_SNR, g = f + ϵ, SNR = 20 * log10(norm(f)/norm(g-f)).

# Input Arguments
- `f::Array{Float64}`: original signal.
- `g::Array{Float64}`: noisy signal.

# Output Argument
- `SNR::Float64`: SNR value.
"""
function compute_SNR(f, g)
    SNR = 20 * log10(norm(f)/norm(g-f))
    return SNR
end

"""
    sort_wavelets(A)

SORT\\_WAVELETS, sort A's column wavelet vectors based on their focused location on Path.

# Input Arguments
- `A::Matrix{Float64}`: whose column vectors are wavelets.

# Output Argument
- `A::Matrix{Float64}`: a matrix with sorted column.
"""
function sort_wavelets(A; order_by_loc = true)
    # sgn = (maximum(A, dims = 1)[:] .> -minimum(A, dims = 1)[:]) .* 2 .- 1
    # A = (A' .* sgn)'

    if order_by_loc
        ord = findmax(abs.(A), dims = 1)[2][:]
        idx = sortperm([j[1] for j in ord])
        A = A[:,idx]
    end

    N = size(A,1)
    sgn = ones(size(A,2))
    mid = Int(round(size(A,2)/2))
    for i in 1:size(A,2)
        cor_res = crosscor(A[:,mid], A[:,i], -Int(ceil(N/2)):Int(floor(N/2)))
        if maximum(cor_res) < -minimum(cor_res)
            sgn[i] = -1
        end
    end
    A = A * Diagonal(sgn)
    return A
end
