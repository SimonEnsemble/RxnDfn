using ProgressMeter
function compute_u(left_bc::Union{Neumann, Dirichlet, Periodic, ConvectiveHeat}, right_bc::Union{Neumann, Dirichlet, Periodic, ConvectiveHeat}, A, discretization::Discretization, st::SpaceTime, u₀, x, sample_time::Float64, f, λ)
    # initialize u at k = 0
    # this u is the value of u(x, t) when t = k - 1 over discretized x points.
    u = u₀.(x)

    # all are unknown in Neumann and ConvectiveHeat, no trimming required
    if isa(left_bc, Dirichlet) && isa(right_bc, Dirichlet)
        u = u[2:end-1]
        x = x[2:end-1]
    elseif isa(left_bc, Periodic)
        u = u[1:end-1]
        x = x[1:end-1]
    else
        if isa(left_bc, Dirichlet)
            u = u[2:end]
            x = x[2:end]
        end
        if isa(right_bc, Dirichlet)
            u = u[1:end-1]
            x = x[1:end-1]
        end
    end

    # initialize u_sample and t_sample to be filled below with specific time steps for graphing purposes
    sample_freq, u_sample, t_sample, sample_counter = sample_variables(st, discretization, sample_time)
    # determine number of unknowns
    nb_unknowns = nb_spatial_unknowns(discretization, left_bc, right_bc)
    # initialize right hand side matrix
    rhs = zeros(Float64, nb_unknowns)
    # take time steps.
    @showprogress 1 "Computing..." for k = 1:discretization.Nₜ
        # time here, inside the loop
        t = discretization.Δt * k

        # Store u and t every sample_time so we can plot.
        if (k == 1) || (k % sample_freq == 0)
            # Different bc's have different num of unknowns, so we account for that here
            # Dirichlet bc's are added later
            if isa(left_bc, Periodic)
                u_sample[1:end-1, sample_counter] = u
            elseif (isa(left_bc, Neumann) || isa(left_bc, ConvectiveHeat)) && (isa(right_bc, Neumann) || isa(right_bc, ConvectiveHeat))
                # accounts for Neumann-Neumann, ConvectiveHeat-ConvectiveHeat, N-CH & CH-N bc's
                u_sample[:, sample_counter] = u
            elseif isa(left_bc, Dirichlet) && isa(right_bc, Dirichlet)
                u_sample[2:end-1, sample_counter] = u
            elseif (isa(left_bc, Neumann) || isa(left_bc, ConvectiveHeat)) && isa(right_bc, Dirichlet)
                u_sample[1:end-1, sample_counter] = u
            elseif isa(left_bc, Dirichlet) && (isa(right_bc, Neumann) || isa(right_bc, ConvectiveHeat))
                u_sample[2:end, sample_counter] = u
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
        if isa(left_bc, Periodic)
            # draw a circle to see.
            # u[-1] is actually u[end] by PBCs
            rhs[1] += λ * (u[end] - 2 * u[1] + u[2])
            # u[end+1] is actually end
            rhs[end] += λ * (u[end - 1] - 2 * u[end] + u[1])
        elseif isa(left_bc, Neumann)
            # -4 and 4 come from making "ghost points" as described by Gustafson (u₀ and uₙ₊₁ are considered our ghost points)
            # ∂ū₀ = (u₂ - u₀)/(2 * Δx) which means u₀ = -2 * Δx * ∂ū₀ + u₂
            # substitute u₀ in the right hand side of the diffusion equation at i = 1 to get the equation below for rhs[1]
            # by the same logic, use uₙ₊₁ = 2 * Δx * ∂ūₗ + uₙ₋₁ for the right hand side of the diffusion equation at i = end
            rhs[1] += λ * (-4.0 * discretization.Δx * left_bc.∂ū + 2 * u[2] - 2 * u[1])
        elseif isa(left_bc, Dirichlet)
            # Simply input the known boundary condition for Dirichlet
            rhs[1] += λ * (left_bc.ū - 2 * u[1] + u[2])
        elseif isa(left_bc, ConvectiveHeat)
            # similar to Neumann, utilizes "ghost" points
            rhs[1] += 2 * λ * (-u[1] * (1 + discretization.Δx * left_bc.K̄) + u[2] + 2 * discretization.Δx * left_bc.K̄ * left_bc.T̄₀)
        end
        if isa(right_bc, Neumann)
            rhs[end] += λ * (4.0 * discretization.Δx * right_bc.∂ū + 2 * u[end - 1] - 2 * u[end])
        elseif isa(right_bc, Dirichlet)
            rhs[end] += λ  * (u[end - 1] - 2 * u[end]  + right_bc.ū)
        elseif isa(right_bc, ConvectiveHeat)
            # similar to Neumann, utilizes "ghost" points
            rhs[end] += 2 * λ * ( u[end] * (discretization.Δx * right_bc.K̄ - 1) + u[end - 1] - 2 * discretization.Δx * right_bc.K̄ * right_bc.T̄₀)
        end

        # solve for u at next time step (overwrite previous)
        u = A \ rhs

    end

    # add in boundary conditions to Dirichlet and Periodic
    if isa(left_bc, Periodic)
        # assign last row to be the same as first since Periodic bc's
        u_sample[end, :] = u_sample[1, :]
        # add last spatial step to x
        push!(x, st.L)
    elseif isa(left_bc, Dirichlet) && isa(right_bc, Dirichlet)
        # insert bc's for first and last rows
        u_sample[1, :] = left_bc.ū
        u_sample[end, :] = right_bc.ū
        # add first and last spatial step to x
        x = vcat([0.0], x, [st.L])
    elseif isa(left_bc, Dirichlet)
        u_sample[1, :] = left_bc.ū
        x=vcat([0.0], x)
    elseif isa(right_bc, Dirichlet)
        u_sample[end, :] = right_bc.ū
        push!(x, st.L)
    end

    return t_sample, x, u_sample
end

compute_u(bc::Periodic, A, discretization::Discretization, st::SpaceTime, u₀, x, sample_time::Float64, f, λ) = compute_u(bc, bc, A, discretization, st, u₀, x, sample_time, f, λ)
