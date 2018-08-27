#Function to create heat map based on x, t, u
import PyPlot; const plt = PyPlot
function draw_heat_map(t, x, u)
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
end
