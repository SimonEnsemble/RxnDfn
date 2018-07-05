#Making array of the exact equation for comparisons
function create_exact_u(t, x, exact_u, D)
    Nₜ = length(t)
    Nₓ = length(x)
    u_e = zeros(Float64, Nₓ, Nₜ)
    for i = 1:Nₓ
        for j = 1:Nₜ
            u_e[i,j] = exact_u(x[i], t[j], D)
        end
    end
    return u_e
end
