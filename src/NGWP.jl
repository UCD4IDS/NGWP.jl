__precompile__()

module NGWP

using LinearAlgebra, SparseArrays, LightGraphs, SimpleWeightedGraphs, Clustering
using JuMP, Clp, Optim, OptimalTransport, MTSG, Statistics
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
export SunFlowerGraph, dualgraph
export pc_ngwp, pairclustering
export vm_ngwp
export ngwp_analysis, ngwp_bestbasis, NGWP_jkl
export scatter_gplot, scatter_gplot!, wiggle, wiggle!
export spike, characteristic, sort_wavelets, transform2D
export getall_expansioncoeffs, approx_error_plot

end
