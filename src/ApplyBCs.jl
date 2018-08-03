#
function apply_bcs!(A, left_bc::Union{Neumann, Dirichlet, ConvectiveHeat}, right_bc::Union{Neumann, Dirichlet, ConvectiveHeat}, λ, discretization)
    if isa(left_bc, Neumann)
        A[1, 2] = -2 * λ
    elseif isa(left_bc, ConvectiveHeat)
        A[1, 1] = 1 + 2 * λ + 2 * λ * discretization.Δx * left_bc.K̄
        A[1, 2] = -2 * λ
    end

    if isa(right_bc, Neumann)
        A[end, end - 1] = -2 * λ
    elseif isa(right_bc, ConvectiveHeat)
        A[end, end] = 1 + 2 * λ - 2 * λ * discretization.Δx * right_bc.K̄
        A[end, end - 1] = -2 * λ
    end

    return A
end
