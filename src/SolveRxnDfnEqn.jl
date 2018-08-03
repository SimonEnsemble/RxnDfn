# solves rxn dfn eqn based off of left and right boundary conditions
function solve_rxn_diffn_eqn(left_bc::Union{Neumann, Dirichlet, ConvectiveHeat}, right_bc::Union{Neumann, Dirichlet, ConvectiveHeat}, f, u₀, D::Float64, Nₓ::Int, st::SpaceTime, sample_time::Float64)
    # ensure that the boundary condition is consistent with the initial condition
    assert_bcs(left_bc, right_bc, u₀, st)
    # discretize time and space, assign λ
    x, λ, discretization = setup_discretization(Nₓ, st, D)
    # determine number of unknowns
    nb_unknowns = nb_spatial_unknowns(discretization, left_bc, right_bc)
    # build tri-diag-matrix without boundary conditions
    A = build_tri_diagonal_matrix(discretization, λ, nb_unknowns)
    # insert boundary condtions into tri-diag-matrix
    apply_bcs!(A, left_bc, right_bc, λ, discretization)
    # computes t, x, and u
    t, x, u = compute_u(left_bc, right_bc, A, discretization, st, u₀, x, sample_time, f, λ)

    return t, x, u
end
