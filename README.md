# JPmesh
Simple Particle Mesh code in Julia for LambdaKCDM cosmology.

It follows the implementation detailed in [K97](https://arxiv.org/abs/astro-ph/9712217) (see also [K17](https://www.worldscientific.com/doi/10.1142/9789813231962_0002)).

## Installation
Install
```julia
using Pkg
Pkg.add(url="https://github.com/tos-1/JPmesh.jl") # for linear P(k) and growth
```

## Example
An example of a run
```julia
using JPmesh
using Plots

# create output directory
output = "data/"
if ! ispath(output) mkdir(output) end
sp = 3 # saves positions every sp steps

# define cosmo params
z  = 30                # initial redshift
ng = 128               # ng^3 particles
nf = 2                 # force resolution: (ng*nf)^3 cells
boxsize = 256.0        # Mpc/h
Om_m = 0.3
Om_l=1-Om_m            # Om_k = 0
simbox = SimBox(Om_m=Om_m, Om_l=Om_l, ng=ng, nf=nf, boxsize=boxsize, z_i=z)
ng = simbox.ng
da = (simbox.a_f-simbox.a_i)/simbox.steps

println("initial redshift z = ", z, " a_i = ", simbox.a_i)
println("number of steps = ", simbox.steps)
println("stepsize = ", da)

# Initial conditions via 2lpt
pos, vel = InitialConditions(simbox);
den = CiC(pos, simbox.ng, simbox.boxsize)
display( heatmap(1:size(den,1), 1:size(den,2), den[:,:,div(ng,2)]) )

# run PM
kickNdrift!(pos, vel, simbox, output, 3)

# plot results
den = CiC(pos, simbox.ng, simbox.boxsize)
den = log.(2.0 .+ den)
display( heatmap( 1:ng, 1:ng, den[:,:,div(ng,2)]) )
```

All the cosmological and simulation parameters are set in the `struct SimBox` defined in JPmesh.jl module. 

![This is an image](https://github.com/tos-1/jpmesh/blob/main/images/LSS.gif)

## License

[MIT](https://choosealicense.com/licenses/mit/)
