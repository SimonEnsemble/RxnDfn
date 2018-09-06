__precompile__(false)

module RxnDfn

include("BoundaryConditions.jl")
include("SpaceTime_Discretization.jl")
include("Assert_BCs.jl")
include("SetupDiscretization.jl")
include("TriDiagMatrix.jl")
include("ApplyBCs.jl")
include("SampleVariables.jl")
include("SpatialUnknowns.jl")
include("Compute_u.jl")
include("SolveRxnDfnEqn.jl")
include("PeriodicSolveRxnDfnEqn.jl")
include("CreateHeatMap.jl")
include("Exact_u.jl")
include("Plot_u.jl")

export BoundaryCondition, Dirichlet, Neumann, Periodic, ConvectiveHeat,
    build_tri_diagonal_matrix, Discretization, SpaceTime, setup_discretization,
    nb_spatial_unknowns, solve_rxn_diffn_eqn, gif_maker,
    draw_heat_map, nb_spatial_unknowns, create_exact_u, apply_bcs!
end
