module JPmesh

import Base.@kwdef
@kwdef struct SimBox
    # cosmo params
    h    :: AbstractFloat = 0.7   
    Om_m :: AbstractFloat = 0.3
    Om_l :: AbstractFloat = 0.7
    Om_b :: AbstractFloat = 0.05
    k_p  :: AbstractFloat = 0.05
    ns :: AbstractFloat = 0.965
    A  :: AbstractFloat = 2.097e-9 # useless
    s8 :: AbstractFloat = 0.80

    # sim parameters
    nf :: Int = 2 # force resolution factor
    ng :: Int = 128
    boxsize :: AbstractFloat = 512 # Mpc/h

    # steps
    z_i :: AbstractFloat = 50  # initial redshift
    a_i :: AbstractFloat = 1/(1+z_i)
    z_f :: AbstractFloat = 0  # final redshift
    a_f :: AbstractFloat = 1/(1+z_f)
    fsteps :: Int = 1.0   # >1, 1 for the minimal number of steps 
    steps :: Int = ceil(Int, fsteps * (a_f-a_i)/(a_i))

    # inital conditions
    seed :: Int = 31415

    # precision of array fields
    prec = Float32  # or Float64 for double 
end
using MatterPower, FFTW, Random, JLD
export SimBox, CiC, trilinear
export D1, EisHu, Fa
export InitialConditions
export kickNdrift!
include("grid.jl")
include("cosmo.jl")
include("initial.jl")
include("pmesh.jl")

end
