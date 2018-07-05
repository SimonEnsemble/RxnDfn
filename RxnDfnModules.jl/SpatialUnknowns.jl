#Determine the number of spatial unknowns at each time step.
function nb_spatial_unknowns(discretization::Discretization, bc::BoundaryCondition)
    if typeof(bc) == Neumann
        return discretization.Nₓ
    elseif typeof(bc) == Dirichlet
        return discretization.Nₓ - 2
    elseif typeof(bc) == Periodic
        return discretization.Nₓ - 1
    end
end
