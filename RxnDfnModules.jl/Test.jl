
push!(LOAD_PATH, "C:\\Users\\Rachel\\Documents\\Research\\RxnDfnModules.jl\\")
using RxnDfn

using Base.Test
@testset "RxnDfn Exact Solution Tests" begin 
    @testset "Dirichlet BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64) 
            return 100 * exp(-t) * sin(π * x)
        end
        D = 1.0
        f(x::Float64, t::Float64, u::Float64) = 100 * (D * π^2 - 1.0) * (e ^ (-t)  * sin(π * x))
        u₀(x::Float64) = 100 * sin(π * x)
        bc = Dirichlet(0.0, 0.0)
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(f, u₀, bc, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
    
    @testset "Neumann BC Test" begin
        # Reaction term
        g(x::Float64) = x^3
        function exact_u(x::Float64, t::Float64, D::Float64) 
            return e^(-π^2 * t) * cos(π * x) + g(x)
        end
        D = 1.0
        f(x::Float64, t::Float64, u::Float64) = -6 * x
        u₀(x::Float64) = cos(π * x) + g(x)
        bc = Neumann(0.0, 3.0)
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(f, u₀, bc, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
    
    @testset "Periodic BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64)
            return exp(-(2 * π)^2 * D * t) * (4 * sin(2 * π * x) + 6 * cos(2 * π * x)) + x * (x - 1)
        end
        D = 0.01
        f(x::Float64, t::Float64, u::Float64) = -2 * D
        u₀(x::Float64) = 4 * sin(2 * π * x) + 6 * cos(2 * π * x) + x * (x - 1)
        bc = Periodic()
        Nₓ = 1000
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1 
        t, x, u = solve_rxn_diffn_eqn(f, u₀, bc, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
end;
