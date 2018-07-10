#Function to create heat map based on x, t, u
function draw_heat_map(t, x, u)
    import PyPlot; const plt = PyPlot
    Nₜ = length(t)
    Nₓ = length(x)

    X = similar(transpose(u))
    T = similar(transpose(u))

    for i = 1:Nₜ
        X[i,:] = x
    end

    for i = 1:Nₓ
        T[:,i] = t
    end

    plt.figure()
    plt.pcolormesh(X,T, transpose(u))
    plt.xlabel("x")
    plt.ylabel("t")
    plt.colorbar(label="Concentration of u")
    plt.savefig("heatmap.png")
    if typeof(bc) == Neumann
        plt.title("Neumann Boundary Conditions")
    elseif typeof(bc) == Dirichlet
        plt.title("Dirichlet Boundary Conditions")
    elseif typeof(bc) == Periodic
        plt.title("Periodic Boundary Conditions")
    end

end
