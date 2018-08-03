function setup_discretization(Nₓ::Int, st::SpaceTime, D::Float64)
    # compute Δx, Nₜ, Δt below.
    discretization = Discretization(Nₓ, NaN, 0, NaN)

    # discretize space.
    x = collect(linspace(0, st.L, discretization.Nₓ)) # includes end points!
    discretization.Δx = x[2] - x[1]
    @printf("%d points in x-discretization. dx = %f\n", discretization.Nₓ, discretization.Δx)

    # discreteize time.
    discretization.Δt = discretization.Δx ^ 2
    discretization.Nₜ = ceil(Int, st.tf / discretization.Δt)
    @printf("%d points in t-discretization. dt = %f\n", discretization.Nₜ, discretization.Δt)

    # nondimensional parameter involved in discreteization
    λ = (D * discretization.Δt) / (2 * (discretization.Δx) ^ 2)

    return x, λ, discretization
end
