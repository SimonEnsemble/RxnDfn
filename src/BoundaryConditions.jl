abstract type BoundaryCondition
end

type Dirichlet <: BoundaryCondition
    ū::Float64
end

type Neumann <: BoundaryCondition
    ∂ū::Float64
end

type Periodic <: BoundaryCondition
end

type ConvectiveHeat <: BoundaryCondition
    T̄₀::Float64
    K̄::Float64 # K̄ = h/k
end
