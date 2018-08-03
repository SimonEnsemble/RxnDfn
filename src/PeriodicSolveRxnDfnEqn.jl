function solve_rxn_diffn_eqn(bc::Periodic, f, u₀, D::Float64, Nₓ::Int, st::SpaceTime, sample_time::Float64)
    # ensure that the boundary condition is consistent with the initial condition
    @assert(u₀(0.0) ≈ u₀(st.L),
        "The initial condition is inconsistent with Periodic boundary conditions.")
    # discretize time and space, assign λ
    x, λ, discretization = setup_discretization(Nₓ, st, D)
    # determine number of unknowns
    nb_unknowns = discretization.Nₓ - 1
    # build tri-diag-matrix without boundary conditions
    A = build_tri_diagonal_matrix(discretization, λ, nb_unknowns)
    # insert boundary condtions into tri-diag-matrix
    A[1, end] = -λ
    A[end, 1] = -λ
    # computes t, x, and u
    t, x, u = compute_u(bc, A, discretization, st, u₀, x, sample_time, f, λ)

    return t, x, u
end
