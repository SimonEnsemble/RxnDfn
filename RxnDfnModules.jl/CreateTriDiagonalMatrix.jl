#Create tridiagonal matrix.

function build_tri_diagonal_matrix(discretization::Discretization, λ::Float64, bc::BoundaryCondition)
    # determine number of unknowns based on boundary condition
    matrix_dim = nb_spatial_unknowns(discretization, bc)

    A = zeros(matrix_dim, matrix_dim)

    for i = 1:matrix_dim
        A[i, i] = 1 + 2 * λ
    end

    for i = 1:matrix_dim - 1
        A[i, i+1] = -λ
        A[i+1, i] = -λ
    end

    if typeof(bc) == Periodic
        A[1, end] = -λ
        A[end, 1] = -λ
    end

    if typeof(bc) == Neumann
        A[1, 2] = -2 * λ
        A[end, end - 1] = -2 * λ
    end

    return A

end
