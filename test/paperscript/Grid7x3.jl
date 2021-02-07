# script for Fig.2, Fig.4, Fig.6

using NGWP, LightGraphs, Plots, LaTeXStrings, MultivariateStats
pyplot(dpi = 200)

Nx, Ny = 7, 3
G = LightGraphs.grid([Nx,Ny]); N = nv(G);
L = Matrix(laplacian_matrix(G))
Q = incidence_matrix(G; oriented = true)
𝛌, 𝚽 = eigen(L); 𝚽 = 𝚽.*sign.(𝚽[1,:])';

# DAG pseudo-metric
distDAG = eigDAG_Distance(𝚽,Q,N)

# MDS embedding into R²
D = distDAG
E = transform(fit(MDS, D, maxoutdim=2, distances=true))

# set up all heatmap plots' positions
dx = 0.01; dy = dx;
xej = zeros(Nx,N); yej=zeros(Ny,N);
a = 5.0; b = 9.0;
for k=1:N
    xej[:,k]=LinRange(E[1,k]-Ny*a*dx,E[1,k]+Ny*a*dx, Nx); yej[:,k]=LinRange(E[2,k]-a*dy,E[2,k]+a*dy, Ny)
end

# find correct 2D index
grid2eig_ind = [1,2,3,6,8,12,15,4,5,7,9,13,16,18,10,11,14,17,19,20,21]; eig2grid_ind = sortperm(grid2eig_ind);
eig2dct = Array{Int64,3}(undef, Nx, Ny, 2); for i in 1:Nx; for j in 1:Ny; eig2dct[i,j,1] = i-1; eig2dct[i,j,2] = j-1; end; end; eig2dct = reshape(eig2dct, (N,2)); eig2dct = eig2dct[eig2grid_ind,:];

#################### Fig.2
plot()
for k=1:N
    heatmap!(xej[:,k],yej[:,k],reshape(𝚽[:,k],(Nx,Ny))',c=:viridis,colorbar=false,ratio=1,annotations=(xej[4,k], yej[3,k]+b*dy, text(latexstring("\\phi_{", string(eig2dct[k,1]), ",", string(eig2dct[k,2]), "}"),10)))
end
plt = plot!(aspect_ratio = 1, xlim = [-1.4, 1.3], ylim = [-1.4, 1.3], grid = false, clims=(-0.4,0.4))
savefig(plt, joinpath(@__DIR__, "../paperfigs/Grid7x3_DAG_MDS.png"))

#################### Fig.4
# first level partition
p1x = [-0.2, 1.0, NaN]; p1y = [1.3, -1.0, NaN]; plot!(p1x, p1y, c = :red, legend = false, width = 3)
# second level partition
p2x = [-1.0, 0.2, NaN, 0.4, 1.2, NaN]; p2y = [-0.8, 0.45, NaN, 0.25, 0.2, NaN]; plot!(p2x, p2y, c=:orange, legend = false, width = 2)
plt = current()
savefig(plt, joinpath(@__DIR__, "../paperfigs/Grid7x3_DAG_2levels_partition.png"))


## Build Dual Graph
W_dual = sparse(dualGraph(distDAG)) #required: sparse dual weighted adjacence matrix

## Assemble NGWPs
ht_elist_PC, ht_vlist_PC = HTree_EVlist(𝚽,W_dual)
wavelet_packet_PC = HTree_wavelet_packet(𝚽,ht_vlist_PC,ht_elist_PC)

ht_elist_VM = ht_elist_PC
wavelet_packet_VM = HTree_wavelet_packet_varimax(𝚽,ht_elist_VM)

## level 2 VM-NGWP vectors
lvl = 3; W_VM = hcat(wavelet_packet_VM[lvl][1], wavelet_packet_VM[lvl][2], wavelet_packet_VM[lvl][3], wavelet_packet_VM[lvl][4]);

wav_kl = [[0 0];[0 1];[0 2];[1 0];[1 1];[1 2];[1 3];[2 0];[2 1];[2 2];[3 0];[3 1];[3 2];[3 3];[2 3];[2 4];[2 5];[2 6];[3 4];[3 5];[3 6]]; wav_kl = wav_kl[eig2grid_ind,:];

reorder_ind = [1,3,2,5,7,4,6,9,10,8,11,13,14,12,16,18,15,17,20,21,19]
W_VM = W_VM[:,reorder_ind[eig2grid_ind]];
sgn = ones(N); sgn[grid2eig_ind[[4,6,8,10,18]]] .= -1; W_VM = W_VM * Diagonal(sgn);

#################### Fig.6
plot()
for k=1:N
    heatmap!(xej[:,k],yej[:,k],reshape(W_VM[:,k],(Nx,Ny))',c=:viridis,colorbar=false,ratio=1,annotations=(xej[4,k], yej[3,k]+b*dy, text(latexstring("\\psi_{", string(wav_kl[k,1]), ",", string(wav_kl[k,2]), "}"),10)))
end
plot!(aspect_ratio = 1, xlim = [-1.4, 1.3], ylim = [-1.4, 1.3], grid = false, clims=(-0.34,0.34))
# first level partition
p1x = [-0.2, 1.0, NaN]; p1y = [1.3, -1.0, NaN]; plot!(p1x, p1y, c = :red, legend = false, width = 3)
# second level partition
p2x = [-1.0, 0.2, NaN, 0.4, 1.2, NaN]; p2y = [-0.8, 0.45, NaN, 0.25, 0.2, NaN]; plot!(p2x, p2y, c=:orange, legend = false, width = 2)
plt = current()
savefig(plt, joinpath(@__DIR__, "../paperfigs/Grid7x3_DAG_VM_NGWP_lvl2_wavelets.png"))
