using Calculus
using ProgressMeter
"""
Arguments:
* f: reaction term f(x, t, u)
* u₀: initial condition u₀(x)
* bc::BoundaryCondition: boundary condition
* D::Float64: diffusion coefficient
* Nₓ::Int: number of spatial discretization points
* st::SpaceTime: space-time over which solution to PDE is approximated
* sample_time::Float64: stores u every sample_time time steps
"""
function solve_rxn_diffn_eqn(f, u₀, bc::BoundaryCondition, D::Float64, Nₓ::Int, st::SpaceTime, sample_time::Float64)

    # compute Δx, Nₜ, Δt below.
    discretization = Discretization(Nₓ, NaN, 0, NaN)

    # ensure that the boundary condition is consistent with the initial condition
    if isa(bc, Dirichlet)
        @assert(bc.ū₀ ≈ u₀(0.0),
            "The initial condition is inconsistent with Dirichlet boundary condition.")
        @assert(isapprox(bc.ūₗ, u₀(st.L), atol = 1e-10),
            "The initial condition is inconsistent with Dirichlet boundary condition.")
    elseif isa(bc, Periodic)
        @assert(u₀(0.0) ≈ u₀(st.L),
            "The initial condition is inconsistent with Periodic boundary conditions.")
    elseif isa(bc, Neumann)
        @assert(isapprox(bc.∂ū₀, derivative(u₀, 0.0), atol = 1e-10),
            "The initial condition is inconsistent with Neumann boundary condition.")
        @assert(isapprox(bc.∂ūₗ, derivative(u₀, st.L), atol = 1e-10),
            "The initial condition is inconsistent with Neumann boundary condition.")
    end

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

    # build tri-diagonal matrix
    A = build_tri_diagonal_matrix(discretization, λ, bc)

    nb_unknowns = nb_spatial_unknowns(discretization, bc)

    # initialize u at k = 0
    # this u is the value of u(x, t) when t = k - 1 over discretized x points.
    u = u₀.(x)

    # trim u and x so that i index corresponds to x and uᵢ,ₖ
    if typeof(bc) == Neumann
        # all are unknown, no trimming required
    elseif typeof(bc) == Dirichlet
        u = u[2:end-1]
        x = x[2:end-1]
    elseif typeof(bc) == Periodic
        u = u[1:end-1]
        x = x[1:end-1]
    end
    @assert(length(u) == nb_unknowns)
    @assert(length(x) == nb_unknowns)

    # initialize right-hand side of matrix eqn. solved at each time step.
    rhs = zeros(Float64, nb_unknowns)

    # initialize u_sample and t_sample to be filled below with specific time steps for graphing purposes
    num_samples = ceil(Int, st.tf / sample_time) # determines size of u_sample and t_sample
    sample_freq = floor(Int, discretization.Nₜ / (num_samples - 1)) # don't count the 0 sample
    u_sample = zeros(Float64, discretization.Nₓ, num_samples) # nb_unknowns = size of u, num_samples, how many u's we will store
    t_sample = zeros(Float64, num_samples)
    sample_counter = 1 # for if loop below to identify column of u_sample, increments below after u_sample is imput

    # take time steps.
    @showprogress 1 "Computing..." for k = 1:discretization.Nₜ
       # time here, inside the loop
        t = discretization.Δt * k

        # Store u and t every so often so we can plot.
        #
        if (k == 1) || (k % sample_freq == 0)
             if typeof(bc) == Neumann
                 u_sample[:, sample_counter] = u
             elseif typeof(bc) == Dirichlet
                u_sample[2:end-1, sample_counter] = u
            elseif typeof(bc) == Periodic
                u_sample[1:end-1, sample_counter] = u
             end
                t_sample[sample_counter] = t - discretization.Δt
                sample_counter += 1
         end

        # build right-hand side vector of matrix eqn.
        for i = 1:nb_unknowns
            # contribution from first order discretization of time derivative
            rhs[i] = u[i]
            # contribution from reaction term
            rhs[i] += f(x[i], t, u[i]) * discretization.Δt
            # TODO contribution from advection?
            # contribution from spatial Laplacian
            if (i != 1) && (i != nb_unknowns)
                rhs[i] += λ * (u[i - 1] - 2 * u[i] + u[i + 1])
            end
         end

         # based on boundary condition, handle first and last component of Laplacian.
         if typeof(bc) == Neumann
         # -4 and 4 come from making "ghost points" as described by Gustafson (u₀ and uₙ₊₁ are considered our ghost points)
         # ∂ū₀ = (u₂ - u₀)/(2 * Δx) which means u₀ = -2 * Δx * ∂ū₀ + u₂
         # substitute u₀ in the right hand side of the diffusion equation at i = 1 to get the equation below for rhs[1]
         # by the same logic, use uₙ₊₁ = 2 * Δx * ∂ūₗ + uₙ₋₁ for the right hand side of the diffusion equation at i = end
             rhs[1] += λ * (-4.0 * discretization.Δx * bc.∂ū₀ + 2 * u[2] - 2 * u[1])
             rhs[end] += λ * (4.0 * discretization.Δx * bc.∂ūₗ + 2 * u[end - 1] - 2 * u[end])
         elseif typeof(bc) == Dirichlet
             # Simply input the known boundary condition for Dirichlet
             rhs[1] += λ * (bc.ū₀ - 2 * u[1] + u[2])
             rhs[end] += λ  * (u[end - 1] - 2 * u[end]  + bc.ūₗ)
         elseif typeof(bc) == Periodic
             # draw a circle to see.
             # u[-1] is actually u[end] by PBCs
             rhs[1] += λ * (u[end] - 2 * u[1] + u[2])
             # u[end+1] is actually end
             rhs[end] += λ * (u[end - 1] - 2 * u[end] + u[1])
         end

         # solve for u at next time step (overwrite previous)
         u = A \ rhs
    end

    if typeof(bc) == Dirichlet
        u_sample[1, :] = bc.ū₀
        u_sample[end, :] = bc.ūₗ
        x = vcat([0.0], x, [st.L])
    elseif typeof(bc) == Periodic
        u_sample[end, :] = u_sample[1, :]
        push!(x, st.L)
    end

    return t_sample, x, u_sample

end
