#Determine the number of spatial unknowns at each time step.
function nb_spatial_unknowns(discretization::Discretization, left_bc::Union{Neumann, Dirichlet, Periodic, ConvectiveHeat}, right_bc::Union{Neumann, Dirichlet, Periodic, ConvectiveHeat})
    if (isa(left_bc, Neumann) || isa(left_bc, ConvectiveHeat)) && (isa(right_bc, Neumann) || isa(right_bc, ConvectiveHeat))
        # accounts for Neumann-Neumann, ConvectiveHeat-ConvectiveHeat, N-CH & CH-N bc's
        return discretization.Nₓ
    elseif isa(left_bc, Dirichlet) && isa(right_bc, Dirichlet)
        return discretization.Nₓ - 2
    else
        # accounts for Periodic, Dirichlet-Neumann, & Dirichlet-ConvectiveHeat bc's
            return discretization.Nₓ - 1
    end
end
