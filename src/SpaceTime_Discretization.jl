mutable struct Discretization
    Nₓ::Int64
    Δx::Float64
    Nₜ::Int64
    Δt::Float64
end

mutable struct SpaceTime
    L::Float64 #spatial extent
    tf::Float64 #time span of simulation
end
