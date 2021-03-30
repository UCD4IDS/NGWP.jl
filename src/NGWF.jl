function nat_spec_filter(l, D; σ = 0.25 * maximum(D), method = :regular, thres = 0.2)
    d = D[:, l]
    𝛍 = exp.(-(d ./ σ).^2)
    𝛍 ./= sum(𝛍)
    if method == :reduced
        𝛍 .*= (𝛍 .> (thres * maximum(𝛍)))
    end
    return 𝛍
end

function ngwf_all_vectors(D; σ = 0.2 * maximum(D))
    𝓤 = zeros(N, 0)
    for l = 1:N
        𝛍 = nat_spec_filter(l, D; σ = σ)
        P = 𝚽 * diagm(𝛍) * 𝚽'
        𝓤 = hcat(𝓤, P)
    end
    return 𝓤
end

function rngwf_all_vectors(D; σ = 0.2 * maximum(D), thres = 0.2)
    𝓤 = zeros(N, 0)
    dic_l2x = Dict()
    for l = 1:N
        𝛍 = nat_spec_filter(l, D; σ = σ, method = :reduced, thres = thres)
        P = 𝚽 * diagm(𝛍) * 𝚽'
        dic_l2x[l] = qr(P, Val(true)).p[1:rank(P, rtol = 10^4 * eps())]
        𝓤 = hcat(𝓤, P[:, dic_l2x[l]])
    end
    return 𝓤, dic_l2x
end

function ngwf_vector(D, l, x; σ = 0.1 * maximum(D))
    P = 𝚽 * diagm(nat_spec_filter(l, D; σ = σ)) * 𝚽'
    ψ = P * spike(x, N)
    ψ ./= norm(ψ, 2)
    return ψ
end

function frame_approx(f, U, V; num_kept = length(f))
    g = U' * f
    ind = sortperm(g.^2; rev = true)[1:num_kept]
    f_approx = zeros(length(f))
    rel_error = [1.0]
    for γ in ind
        f_approx += g[γ] * V[:, γ]
        f_res = f - f_approx
        push!(rel_error, norm(f_res)/norm(f))
    end
    return rel_error, f_approx
end


function rngwf_lx(dic_l2x)
    N = length(dic_l2x)
    Γ = (Tuple{Int64,Int64})[]
    for l = 1:N
        for x in dic_l2x[l]
            push!(Γ, (l, x))
        end
    end
    return Γ
end
