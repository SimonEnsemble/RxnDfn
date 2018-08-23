# ReactionDiffusionEqn
![Diffusion Equation](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/DfnEqnPic.PNG)

ReactionDiffusionEqn is a Julia package that can solve a reaction diffusion equation based on the following input:
1. Diffusion Coefficient (D)
2. Reaction Term (f(x, t, u))
3. Initial Condition (u₀)
4. Boundary Conditions (Dirichlet, Neumann, Convective Heat (aka Robin), or Periodic)
  ```Julia
  bc = Dirichlet(ū::Float64) # boundary condition must be specified
  bc = Neumann(∂ū::Float64) # boundary condition derivative must be specified
  bc = Periodic() # nothing needs to be specified
  bc = ConvectiveHeat(T̄₀::Float64, K̄::Float64) # ambient temperature (T̄₀) and thermal conductivity (K̄) must be specified
  ```
5. Number of Spatial Steps (Nₓ - number of spatial discretization points)
6. Space Time (space-time over which solution to PDE is approximated)
  ```Julia
  st = SpaceTime(L::Float64, tf::Float64)  # L is the spatial extent and tf is the time span of simulation)
  ```
7. Sample Time (stores u every sample_time time steps)

We approximate the solution to the PDE numerically using a finite difference spatial and temporal discretization with the Crank-Nicolson Method:
![Crank-Nicolson Equation](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/Crank-NicolsonEqnPic.PNG)

Approximating this complex equation is as easy as passing eight (seven if periodic, see periodic example for further explanation) variables into the solve_rnx_diffn_eqn function.
Below is an example with Dirichlet Boundary Conditions and the corresponding heat map:
```Julia
D = 1.0
f(x::Float64, t::Float64, u::Float64) = 100 * (D * π^2 - 1.0) * (e ^ (-t)  * sin(π * x))
u₀(x::Float64) = 100 * sin(π * x)
left_bc = Dirichlet(0.0)
right_bc = Dirichlet(0.0)
Nₓ = 100
st = SpaceTime(1.0, 1.0)
sample_time = 0.1

t, x, u = solve_rxn_diffn_eqn(left_bc, right_bc, f, u₀, D, Nₓ, st, sample_time)

draw_heat_map(t, x, u)
```
![Dirichlet Heat Map](https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/DirichletHeatMap.png)

ReactionDiffusionEqn can also handle mixed boundary conditions, such as Dirichlet-Neumann, Dirichlet-ConvectiveHeat, Neumann-ConvectiveHeat, etc.

The following functions are also available to aid in visualizations:
```Julia
draw_heat_map(t, x, u)
```
As shown in the example above, the draw_heat_map function outputs a heat map using PyPlot.

```Julia
gif_maker(t, x, u, st, sample_time)
```
The gif_maker function plots the approximated u at each timestep.

<img src="https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/DirichletGif.gif" width="600" height="400" title="DirichletGif">

## Examples

### Dirichlet Boundary Conditions
<img src="https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/DirichletHeatTransferPic.png" width="440" height="329" title="Dirichlet Heat Transfer">

Use the ReactionDiffusionEqn package to determine how the temperature will change over time in an insulated, one-dimensional rod when it is heated at a particular location while the ends of the rod are held at a fixed temperature.

### Neumann Boundary Conditions
<img src="https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/NeumannLakePic.png" width="473" height="237" title="Neumann Fish Boundary Conditions">

Imagine you have a one dimensional lake with a town at one end and a forest at the other as pictured above. We can use the following equation to model the fish population in the lake if people fish there:   

<img src="https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/FishHarvestingEqn.PNG" width="357" height="96" title="Fish Harvesting Reaction Diffusion Equation">

The reaction term consists of the logistic growth model to simulate the population density of the fish and `ch(x)p(u)` to model the harvesting of the fish.

:fish: `c` is the harvesting rate.

:fish: `h(x)` is the harvesting distribution (i.e. how far from the town the person fishing drops a line in the lake).

:fish: `p(u)` is the probability of catching a fish based on the population density of the fish.

### Periodic Boundary Conditions
<img src="https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/GrassShadeExample__NoBackground.png" width="730" height="332" title="Periodic Grassland Model Example">

One can model the changes in grass density where a tree shades some grass and fauna eat the grass with the following equation:

<img src="https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/GrassEquation.PNG" width="357" height="96" title="Grass Growth Reaction Diffusion Equation">

The reaction term consists of the logistic growth model to simulate the density of the grass and `- mb - hqb` to model death of the grass.

:elephant: `g` is the growth rate of the grass

:elephant: `m` is the mortality rate of the grass

:elephant: `h` is the biomass density of the fauna

:elephant: `q` is the rate of consumption of the grass

When using Periodic boundary conditions, there is no need for a left and right boundary condition. Only one boundary condition is needed:
```Julia
bc = Periodic()
t, x, u = solve_rxn_diffn_eqn(bc, f, u₀, D, Nₓ, st, sample_time)
```

### Dirichlet-Convective Heat Boundary Conditions
<img src="https://github.com/SimonEnsemble/RxnDfn/blob/master/Images/DirichletConvective.png" width="593" height="322" title="Dirichlet-Convective Heat Example">

The ReactionDiffusionEqn software can also handle mixed boundary conditions. In this example, we mix Dirichlet and Convective Heat boundary conditions. We have an insulated, one-dimensional rod and can evaluate the change in temperature as the rod is heated while one end of the rod is held at a fixed temperature.
