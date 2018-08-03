#Create tridiagonal matrix.
function build_tri_diagonal_matrix(discretization::Discretization, λ::Float64, nb_unknowns::Int)
    A = zeros(nb_unknowns, nb_unknowns)

    for i = 1:nb_unknowns
        A[i, i] = 1 + 2 * λ
    end

    for i = 1:nb_unknowns - 1
        A[i, i+1] = -λ
        A[i+1, i] = -λ
    end

    return A
end
