"""
    best_basis_selection(f, wavelet_packet, parent_pointer)

BEST BASIS SELECTION inorder to approximate graph signal f, it selects best basis from wavelet\\_packet based on ℓ₁-norm of coefficients.

# Input Arguments
- `f::Array{Float64}`: graph signal.
- `wavelet_packet::Array{Array{Matrix{Float64}}}`: a wavelet packet tree.
- `parent_pointer::Array{Array{Array{Int64}}}`: the parent pointer for the tree.

# Output Argument
- `Wav::Matrix{Float64}`: columns are a set of wavelet ONB.

"""
function best_basis_selection(f, wavelet_packet, parent_pointer)
    ht_coeff_L1 = HTree_coeff_wavelet_packet(f,wavelet_packet)[2]
    bb_locs = best_basis_algorithm(ht_coeff_L1, parent_pointer)
    Wav = assemble_wavelet_basis(bb_locs, wavelet_packet)
    return Wav
end

function best_basis_algorithm(ht_coeff_L1, parent)
    Lvl = length(parent)
    dvec = [[Lvl+1, k] for k = 1:length(ht_coeff_L1[end])]
    Subspaces = [[_ for _ in pair] for pair in parent[Lvl]]

    for lvl = Lvl:-1:1
        dvec_copy = copy(dvec)
        Subspaces_copy = copy(Subspaces)
        tmp_count = 0
        tmp_loc = []
        for i in 1:length(Subspaces)
            pair = Subspaces[i]
            if length(pair) == 1
                ind = findall(x -> x == dvec[pair[1]],dvec_copy)[1]
                deleteat!(dvec_copy,ind)
                insert!(dvec_copy,ind, [lvl, i])
                Subspaces_copy[i] = [ind]
            else
                if compute_subspace_cost(dvec,pair,ht_coeff_L1) > ht_coeff_L1[lvl][i]
                    ind = findall(x -> x == dvec[pair[1]],dvec_copy)[1]
                    deleteat!(dvec_copy,ind:ind+length(pair)-1)
                    insert!(dvec_copy,ind,[lvl, i])

                    tmp_count -= length(pair) - 1
                    push!(tmp_loc, i)

                    Subspaces_copy[i] = [ind]
                else
                    Subspaces_copy[i] = [k + tmp_count for k in Subspaces[i]]
                end
            end
        end

        dvec = dvec_copy
        if lvl > 1
            Subspaces = [union_array_of_arrays(Subspaces_copy[pair]) for pair in parent[lvl-1]]
        else
            Subspaces = [[1]]
        end
    end

    return dvec
end

"""
    best_basis_selected_subspaces(f, wavelet_packet, ht_elist)
Find best basis selected subspaces V^{⋆ (j)}_k. The resulting subspaces, i.e. (j+1,k+1)s, are ordered from the left of the binary tree to the right.
"""
function best_basis_selected_subspaces(f, wavelet_packet, ht_elist)
    parent_pointer = HTree_findParent(ht_elist)
    ht_coeff_L1 = HTree_coeff_wavelet_packet(f,wavelet_packet)[2]
    subspace_locs = best_basis_algorithm(ht_coeff_L1, parent_pointer)
    return subspace_locs
end

"""
    NGWP_jkl(f, wavelet_packet, ht_elist)
Find the best basis corresponding (j,k,l).
"""
function NGWP_jkl(f, wavelet_packet, ht_elist)
    bb_locs = best_basis_selected_subspaces(f, wavelet_packet, ht_elist)
    ind2jkl = (Tuple{Int64,Int64,Int64})[]
    for i in 1:length(bb_locs)
        j = bb_locs[i][1]-1
        k = bb_locs[i][2]-1
        for l in 1:length(ht_elist[j][k+1])
            push!(ind2jkl, (j,k,l-1))
        end
    end
    return ind2jkl
end

function compute_subspace_cost(dvec,arr,ht_coeff_L1)
    s = 0
    for ele in dvec[arr]
        s += ht_coeff_L1[ele[1]][ele[2]]
    end
    return s
end

function union_array_of_arrays(arr)
    s = []
    for k in arr
        s = union(s, k)
    end
    return s
end
