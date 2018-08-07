using RxnDfn

using Base.Test
@testset "RxnDfn Exact Solution Tests" begin
     @testset "Periodic BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64)
            return exp(-(2 * π)^2 * D * t) * (4 * sin(2 * π * x) + 6 * cos(2 * π * x)) + x * (x - 1)
        end
        D = 0.01
        f(x::Float64, t::Float64, u::Float64) = -2 * D
        u₀(x::Float64) = 4 * sin(2 * π * x) + 6 * cos(2 * π * x) + x * (x - 1)
        bc = Periodic()
        Nₓ = 100
        st = SpaceTime(1.0, 0.1)
        sample_time = 0.01
        t, x, u = solve_rxn_diffn_eqn(bc, f, u₀, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;

    @testset "Dirichlet BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64)
            return 100 * exp(-t) * sin(π * x)
        end
        D = 1.0
        f(x::Float64, t::Float64, u::Float64) = 100 * (D * π^2 - 1.0) * (e ^ (-t)  * sin(π * x))
        u₀(x::Float64) = 100 * sin(π * x)
        left_bc = Dirichlet(0.0)
        right_bc = Dirichlet(0.0)
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(left_bc, right_bc, f, u₀, D, Nₓ, st, sample_time)
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
        left_bc = Neumann(0.0)
        right_bc = Neumann(3.0)
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(left_bc, right_bc, f, u₀, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
#=    #TODO ConvectiveHeat BCs
    @testset "ConvectiveHeat BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64)
            return
        end
        D = 1.0
        f(x::Float64, t::Float64, u::Float64) =
        u₀(x::Float64) =
        left_bc = ConvectiveHeat( , )
        right_bc = ConvectiveHeat( , )
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(left_bc, right_bc, f, u₀, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
 =#  
    @testset "Dirichlet-Neumann BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64)
            return e^(-(π^2 * t) / 4) * sin((π * x) / 2)
        end
        D = 1.0
        f(x::Float64, t::Float64, u::Float64) = 0.0
        u₀(x::Float64) = sin((π * x) / 2)
        left_bc = Dirichlet(0.0)
        right_bc = Neumann(0.0)
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(left_bc, right_bc, f, u₀, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
#=    TODO Dirichlet-ConvectiveHeat BC
    @testset "Dirichlet-ConvectiveHeat BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64)
            return
        end
        D = 1.0
        f(x::Float64, t::Float64, u::Float64) =
        u₀(x::Float64) =
        left_bc = Dirichlet()
        right_bc = ConvectiveHeat( , )
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(left_bc, right_bc, f, u₀, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
       TODO Neumann-ConvectiveHeat BC
    @testset "Neumann-ConvectiveHeat BC Test" begin
        function exact_u(x::Float64, t::Float64, D::Float64)
            return
        end
        D = 1.0
        f(x::Float64, t::Float64, u::Float64) =
        u₀(x::Float64) =
        left_bc = Neumann()
        right_bc = ConvectiveHeat( , )
        Nₓ = 100
        st = SpaceTime(1.0, 1.0)
        sample_time = 0.1
        t, x, u = solve_rxn_diffn_eqn(left_bc, right_bc, f, u₀, D, Nₓ, st, sample_time)
        u_e = create_exact_u(t, x, exact_u, D)
        @test isapprox(u, u_e, rtol=.001)
    end;
=#
end;
