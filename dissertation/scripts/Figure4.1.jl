cd(@__DIR__); include("setups/distROT_distTSD_ratio.jl")
gr(dpi = 200)

## (a) P64
include("setups/path64.jl")
ρ1 = generate_ROT_TSD_ratio(nsim, 𝚽, ∇𝚽, 𝛌, Q)
plt = ROT_TSD_ratio_histogram(ρ1)
savefig(plt, "../figs/Path64_ROT_TSD.png")

## (b) P₇ x P₃
include("setups/grid7x3.jl")
ρ2 = generate_ROT_TSD_ratio(nsim, 𝚽, ∇𝚽, 𝛌, Q)
plt = ROT_TSD_ratio_histogram(ρ2)
savefig(plt, "../figs/Grid7x3_ROT_TSD.png")

## (c) Erdos Rényi
include("setups/er.jl")
ρ3 = generate_ROT_TSD_ratio(nsim, 𝚽, ∇𝚽, 𝛌, Q)
plt = ROT_TSD_ratio_histogram(ρ3)
savefig(plt, "../figs/Erdos_Renyi_ROT_TSD.png")

## (d) weighted RGC100
include("setups/rgc100.jl")
ρ4 = generate_ROT_TSD_ratio(nsim, 𝚽, ∇𝚽, 𝛌, Q; edge_length = edge_length)
plt = ROT_TSD_ratio_histogram(ρ4)
savefig(plt, "../figs/wRGC100_ROT_TSD.png")

## Table 4.1
display_basic_stats([ρ1, ρ2, ρ3, ρ4])
