abstract type BoundaryCondition
end

mutable struct Dirichlet <: BoundaryCondition
    ū::Float64
end

mutable struct Neumann <: BoundaryCondition
    ∂ū::Float64
end

mutable struct Periodic <: BoundaryCondition
end

mutable struct ConvectiveHeat <: BoundaryCondition
    T̄₀::Float64
    K̄::Float64 # K̄ = h/k
end
