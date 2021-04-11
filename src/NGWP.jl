__precompile__()

module NGWP

using LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs, Clustering
using JuMP, Clp, Optim, OptimalTransport, MTSG, Statistics, QuadGK, Arpack
import Plots: plot, plot!, scatter, scatter!
import StatsBase: crosscor

filenames = readdir(@__DIR__)

for f in filenames
    if f == "NGWP.jl"
        continue
    elseif occursin(r"^.*\.jl$",f)
        include(joinpath(@__DIR__, f))
    end
end

export eigDAG_Distance, eigHAD_Distance, eigHAD_Affinity
export eigROT_Distance, ROT_Distance, eigEMD_Distance
export eigsROT_Distance
export eigTSD_Distance, K_functional
export natural_eigdist
export SunFlowerGraph, dualgraph
export pc_ngwp, pairclustering
export vm_ngwp, varimax
export lp_ngwp, rising_cutoff, find_pairinds, pair_inds_shadding
export LPHGLET_Synthesis, LPHGLET_Analysis_All, HGLET_dictionary, LPHGLET_dictionary
export unitary_folding_operator, keep_folding!
export ngwp_analysis, ngwp_bestbasis, NGWP_jkl
export natural_eigdist
export nat_spec_filter, ngwf_all_vectors, rngwf_all_vectors, ngwf_vector, frame_approx, rngwf_lx
export scatter_gplot, scatter_gplot!, wiggle, wiggle!
export standardize_eigenvectors!, spike, characteristic, χ, sort_wavelets, transform2D
export getall_expansioncoeffs, approx_error_plot, getall_expansioncoeffs2, approx_error_plot2

end
