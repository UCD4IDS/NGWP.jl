__precompile__()

module NGWP

using LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs, Clustering, JuMP, Clp, Optim, OptimalTransport, MTSG, Statistics
import Plots: plot, plot!, scatter, scatter!
import StatsBase:crosscor

filenames = readdir(@__DIR__)

for f in filenames
    if f == "NGWP.jl"
        continue
    elseif occursin(r"^.*\.jl$",f)
        include(joinpath(@__DIR__, f))
    end
end

export SunFlowerGraph, dualGraph, HTree_EVlist, HTree_wavelet_packet, HTree_wavelet_packet_varimax, scatter_gplot, scatter_gplot!, signal_transform_coeff, approx_error_plot2, HTree_findParent, best_basis_selection, eigDAG_Distance, eigHAD_Distance, eigROT_Distance, eigTSD_Distance, spike, characteristic, sort_wavelets, wiggle, wiggle!

end
