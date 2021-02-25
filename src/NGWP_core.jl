
"""
    ngwp_analysis(G::GraphSig, wavelet_packet::Array{Float64,3})

For a GraphSig object `G`, generate the matrix of NGWP expansion coefficients.

# Input Arguments
- `G::GraphSig`: an input GraphSig object
- `wavelet_packet::Array{Float64,3}`: the varimax wavelets packet.

# Output Argument
- `dmatrix::Array{Float64,3}`: the expansion coefficients matrix.

"""
function ngwp_analysis(G::GraphSig, wavelet_packet::Array{Float64,3})
    f = G.f
    fcols = size(f, 2)
    (N, jmax, ) = Base.size(wavelet_packet)

    dmatrix = zeros(N, jmax, fcols)
    dmatrix[:, 1, :] = f

    for j = 2:jmax
        for i = 1:N
            dmatrix[i, j, :] = f' * wavelet_packet[i, j, :]
        end
    end

    return dmatrix
end


"""
    const_proj_wavelets(𝚽,vlist,elist; method = "Modified Gram-Schmidt with Lp Pivoting")

construct projection wavelets, i.e., project δ ∈ vlist onto span({φⱼ| j ∈ elist}).

# Input Arguments
- `𝚽::Matrix{Float64}`: graph Laplacian eigenvectors 𝚽
- `vlist::Array{Int}`: the list of considered node indices.
- `elist::Array{Int}`: the list of considered eigenvector indices.
- `method::Symbol`: default is `:MGSLp`. other options: `:IP` (Iterative-Projection),
    `GS` (Gram Schmidt).

# Output Argument
- `Wav::Matrix{Float64}`: a matrix whose columns are projected wavelet vectors.

"""
function const_proj_wavelets(𝚽, vlist, elist; method = :MGSLp)
    if length(vlist) == 1
        return 𝚽[:, elist]
    end
    N = size(𝚽, 1)
    m = length(vlist)
    Wav = zeros(N, m)

    B = 𝚽[:, elist]

    if method == :IP
        for k in 1:length(vlist)
            wavelet = Proj(spike(vlist[k], N), B)
            Wav[:, k] .= wavelet ./ norm(wavelet)
            B = wavelet_perp_Matrix(wavelet, B)
        end
    elseif method == :GS
        P = B * B'
        for k in 1:length(vlist)
            wavelet = P * spike(vlist[k], N)
            Wav[:,k] .= wavelet ./ norm(wavelet)
        end
        Wav, complement_dim = gram_schmidt(Wav)
        if complement_dim != 0
            complement_space = B * nullspace(Wav' * B)
            Wav = hcat(Wav, complement_space)
        end
    elseif method == :MGSLp
        P = B * B'
        for k in 1:length(vlist)
            wavelet = P * spike(vlist[k], N)
            Wav[:,k] .= wavelet ./ norm(wavelet)
        end
        Wav, complement_dim = mgslp(Wav)
        if complement_dim != 0
            complement_space = B * nullspace(Wav' * B)
            Wav = hcat(Wav, complement_space)
        end
    else
        error("Do not support method ", method)
    end

    return Wav
end

"""
    function NGWP_jkl(GP_star::GraphPart, drow::Int, dcol::Int)

Generate the (j,k,l) indices for the NGWP basis vector corresponding to the coefficient dmatrix(drow,dcol)

### Input Arguments
* `GP_star::GraphPart`: a GraphPart object of the dual grpah
* `drow::Int`: the row of the expansion coefficient
* `dcol::Int`: the column of the expansion coefficient

### Output Argument
* `j`: the level index of the expansion coefficient
* `k`: the subregion in dual graph's index of the expansion coefficient
* `l`: the tag of the expansion coefficient
"""
function NGWP_jkl(GP_star::GraphPart, drow::Int, dcol::Int)
    (j, k, l) = GHWT_jkl(GP_star, drow, dcol)
    return j,k,l
end
