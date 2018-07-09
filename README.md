# ReactionDiffusionEqn
![Diffusion Equation](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/DfnEqnPic.PNG)

ReactionDiffusionEqn is a Julia package that can solve a reaction diffusion equation based on the following input:
1. Reaction Term (f(x, t, u))
2. Initial Condition (u₀)
3. Boundary Conditions (Dirichlet, Neumann, or Periodic)
4. Diffusion Coefficient (D)
5. Number of Spatial Steps (Nₓ - number of spatial discretization points)
6. Space Time (space-time over which solution to PDE is approximated)
7. Sample Time (stores u every sample_time time steps)

We approximate the solution to the PDE numerically using a finite difference spatial and temporal discretization with the Crank-Nicolson Method:
![Crank-Nicolson Equation](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/Crank-NicolsonEqnPic.PNG)

Approximating this complex equation is as easy as passing in the seven variables into the solve_rnx_diffn_eqn function:
```Julia
f(x::Float64, t::Float64, u::Float64) = 100 * (D * π^2 - 1.0) * (e ^ (-t)  * sin(π * x))
u₀(x::Float64) = 100 * sin(π * x)
bc = Dirichlet(0.0, 0.0)
D = 1.0
Nₓ = 20
st = SpaceTime(1.0, 1.0)
sample_time = 0.1

t, x, u = solve_rxn_diffn_eqn(f, u₀, bc, D, Nₓ, st, sample_time)
```

The following functions are also available to aid in visualizations:
```Julia
gif_maker(t, x, u)
```
The gifmaker function plots the approximated u at each timestep.

```Julia
draw_heat_map(t, x, u)
```
The draw_heat_map function outputs a heat map using PyPlot.

### Dirichlet Boundary Conditions
![Dirichlet Heat Transfer](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/DirichletHeatTransferPic.png)

### Neumann Boundary Conditions
![Neumann Fish Boundary Conditions](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/NeumannFishPic.png)

### Periodic Boundary Conditions
![Periodic Thin Ring Example](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/PeriodicRingPic.png)
