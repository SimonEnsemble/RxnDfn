module RxnDfn

include("BoundaryConditions.jl")
include("SpaceTime_Discretization.jl")
include("SpatialUnknowns.jl")
include("CreateTriDiagonalMatrix.jl")
include("SolveRxnDfnEqn.jl")
include("Plot_u.jl")
include("CreateHeatMap.jl")
include("Exact_u.jl")

export BoundaryCondition, Dirichlet, Neumann, Periodic,
    build_tri_diagonal_matrix, Discretization, SpaceTime,
    nb_spatial_unknowns, solve_rxn_diffn_eqn, gif_maker,
    draw_heat_map, create_exact_u

end
