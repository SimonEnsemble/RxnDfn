# ensure that the boundary condition is consistent with the initial condition
using Calculus
function assert_bcs(left_bc::Union{Neumann, Dirichlet, ConvectiveHeat}, right_bc::Union{Neumann, Dirichlet, ConvectiveHeat}, u₀, st::SpaceTime)
    if isa(left_bc, Dirichlet)
        @assert(left_bc.ū ≈ u₀(0.0),
            "The initial condition is inconsistent with Dirichlet boundary condition.")
    elseif isa(left_bc, Neumann)
        @assert(isapprox(left_bc.∂ū, derivative(u₀, 0.0), atol = 1e-10),
            "The initial condition is inconsistent with Neumann boundary condition.")
    elseif isa(left_bc, ConvectiveHeat)
        @assert(isapprox(K̄ * u₀(0.0) - T̄₀, derivative(u₀, 0.0), atol = 1e-10))
    end

    if isa(right_bc, Dirichlet)
        @assert(isapprox(right_bc.ū, u₀(st.L), atol = 1e-10),
            "The initial condition is inconsistent with Dirichlet boundary condition.")
    elseif isa(right_bc, Neumann)
        @assert(isapprox(right_bc.∂ū, derivative(u₀, st.L), atol = 1e-10),
            "The initial condition is inconsistent with Neumann boundary condition.")
    elseif isa(right_bc, ConvectiveHeat)
        @assert(isapprox(K̄ * u₀(st.L) - T̄₀, derivative(u₀, st.L), atol = 1e-10))
    end
end
