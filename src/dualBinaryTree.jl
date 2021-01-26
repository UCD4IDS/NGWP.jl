"""
    dualBinaryTree(ptr, N)
Draw the dual domain's binary partition tree.
"""
function dualBinaryTree(ptr, N)
    n = 0 # number of nodes in the tree
    jmax = length(ptr)+1
    for j = 1:jmax-1
        n += length(ptr[j])
    end
    n += N

    node_ind = (Array{Int64,1})[]
    ind = 1
    for j = 1:jmax-1
        lvl_nodes = (Int64)[]
        for i = 1:length(ptr[j])
            push!(lvl_nodes, ind)
            ind += 1
        end
        push!(node_ind, lvl_nodes)
    end
    push!(node_ind, [ind-1+i for i=1:N])

    BT = Graph(n)
    X = zeros(n,2)
    hv = ["(0,0)" for i = 1:n]
    for j = 1:jmax-1
        for i in 1:length(ptr[j])
            c1 = ptr[j][i][1]
            add_edge!(BT, Edge(node_ind[j][i], node_ind[j+1][c1]))
            X[node_ind[j+1][c1],:] = X[node_ind[j][i],:] + [0,1]
            hv[node_ind[j+1][c1]] = "($(j),$(c1-1))"
            if length(ptr[j][i]) == 2
                c2 = ptr[j][i][2]
                add_edge!(BT, Edge(node_ind[j][i], node_ind[j+1][c2]))
                X[node_ind[j+1][c1],:] = X[node_ind[j][i],:] + [-2^(jmax-j),1]
                X[node_ind[j+1][c2],:] = X[node_ind[j][i],:] + [2^(jmax-j),1]
                hv[node_ind[j+1][c2]] = "($(j),$(c2-1))"
            end
        end
    end
    return X, BT, node_ind, hv
end
