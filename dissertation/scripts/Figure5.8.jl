cd(@__DIR__); include("setups/simpletree.jl");
pyplot(dpi = 200)

## (a)
plot(size = (200,500))
gplot!(W, X, width = 1)
scatter_gplot!(X[ib1, :]; c = :pink)
scatter_gplot!(X[ib2, :]; c = :orange)
scatter_gplot!(X[ib3, :]; c = :green)
scatter_gplot!(X[ib4, :]; c = :yellow)
scatter_gplot!(X[ijc, :]; c = :red)
scatter_gplot!(X[ir, :]; c = :grey)
plot!(framestyle = :none)
plot!([0, 0], [0.3, 2], line = :arrow, c = :black, lw = 1)
plt = annotate!(0, 3, text("root", 9))
savefig(plt, "../figs/SimpleTree.png")


## (b)
plot(size = (200, 500))
gplot!(Ws, Xs, width = 1)
scatter_gplot!(Xs[11:11, :]; c = :pink, ms = 4)
scatter_gplot!(Xs[10:10, :]; c = :orange, ms = 4)
scatter_gplot!(Xs[13:13, :]; c = :green, ms = 4)
scatter_gplot!(Xs[12:12, :]; c = :yellow, ms = 4)
scatter_gplot!(Xs[[2,4,6,8], :]; c = :red, ms = 4)
scatter_gplot!(Xs[[1,3,5,7,9], :]; c = :grey, ms = 4)
plt = plot!(framestyle = :none, xlims = [minimum(X[:, 1]), maximum(X[:, 1])],
            ylims = [minimum(X[:, 2]), maximum(X[:, 2])])
savefig(plt, "../figs/SimpleTree_simplified.png")
