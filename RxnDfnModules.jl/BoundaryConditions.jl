abstract type BoundaryCondition
end

type Dirichlet <: BoundaryCondition
    ū₀::Float64
    ūₗ::Float64
end

type Neumann <: BoundaryCondition
    ∂ū₀::Float64
    ∂ūₗ::Float64
end

type Periodic <: BoundaryCondition
end
