#################################################################
# particle mesh code. 
# References:
# K97: https://arxiv.org/abs/astro-ph/9712217
# K17: https://www.worldscientific.com/doi/10.1142/9789813231962_0002 
#################################################################

"solve for poisson potential: computes the reduced gravitational potential starting
 from the particles' distribution."
function poisson(pos::Vector{Vector{T}}, a::AbstractFloat, simbox::SimBox)::Array{T,3} where{T<:AbstractFloat}
    ng = simbox.ng*simbox.nf
    boxsize = simbox.boxsize
    # RHS of Eq.(9) of K97
    den = CiC(pos, ng, boxsize ) * simbox.Om_m * 1.5 / a
    dk = rfft(den)

    ng02 = div(ng,2)

    # Eq.(25) of K17, dk the potential
    for i in 1:ng02+1, j in 1:ng, k in 1:ng
       kx = if i<ng02+1 i-1 else (i-1-ng) end
       ky = if j<ng02+1 j-1 else (j-1-ng) end
       kz = if k<ng02+1 k-1 else (k-1-ng) end
       Gk = 0.5 * (boxsize/ng)^2  * ( cos(pi * kx/ng02) + cos(pi * ky/ng02) + cos(pi * kz/ng02) -3 )^(-1)
       dk[i,j,k] *= Gk
    end
    dk[1,1,1] = 0.0
    irfft(dk, ng)
end


"find acceleration from gravitational potential using finite difference approximation 2nd order.
 Interpolate the acceleration on the particles' positions using trilinear."
function newton(phi::Array{T,3}, pos::Vector{Vector{T}}, simbox::SimBox)::Vector{Vector{T}} where{T<:AbstractFloat}
    ng = size(phi)[1] 
    boxsize = simbox.boxsize
    cell = boxsize/ng
    Fx = similar(phi)
    Fy = similar(phi)
    Fz = similar(phi)
    
    for i in 1:ng, j in 1:ng, k in 1:ng
      Fx[i,j,k] = -(2/3 * (phi[1+mod(i-1+1,ng),j,k] - phi[1+mod(i-1-1,ng),j,k]) - 1/12 * (phi[1+mod(i-1+2,ng),j,k] - phi[1+mod(i-1-2,ng),j,k])) / cell
      Fy[i,j,k] = -(2/3 * (phi[i,1+mod(j-1+1,ng),k] - phi[i,1+mod(j-1-1,ng),k]) - 1/12 * (phi[i,1+mod(j-1+2,ng),k] - phi[i,1+mod(j-1-2,ng),k])) / cell
      Fz[i,j,k] = -(2/3 * (phi[i,j,1+mod(k-1+1,ng)] - phi[i,j,1+mod(k-1-1,ng)]) - 1/12 * (phi[i,j,1+mod(k-1+2,ng)] - phi[i,j,1+mod(k-1-2,ng)])) / cell
    end
    # assign acceleration to particles
    ax = trilinear(Fx, pos, boxsize)
    ay = trilinear(Fy, pos, boxsize)
    az = trilinear(Fz, pos, boxsize)

    return [ax, ay, az]
end

"integrate using kick-drift-kick scheme, save pos to output folder, e.g. 'data/'"
function kickNdrift!(pos::Vector{Vector{T}},vel::Vector{Vector{T}}, simbox::SimBox, output::String, sp::Int=1) where{T<:AbstractFloat}

    boxsize = convert(simbox.prec, simbox.boxsize)
    Om_m = simbox.Om_m
    Om_l = simbox.Om_l
    h = simbox.h
    da = (simbox.a_f-simbox.a_i)/simbox.steps
    a  = simbox.a_i
    save(output * "pos_0.jld", "pos", pos, "z", round(simbox.z_i,digits=2))
    vel .*= (simbox.a_i-da/2) # rewrite as momenta, eq(6,8) K97, x0H0=1
    for s in 1:simbox.steps
        println("step= ",s)
        phi = poisson(pos, a, simbox)
        acc = newton(phi, pos, simbox)
        vel .+= acc .* Fa(a,Om_l,Om_m) * da        # kick
        a += da/2
        pos .+= vel .* Fa(a,Om_l,Om_m) / a^2 * da  # drift
	a += da/2
	PBC!(pos,boxsize)
	if s%sp==0 save(output * "pos_" * string(s) * ".jld", "pos", pos, "z", round(1/a-1,digits=2)) end
    end
    phi = poisson(pos, a, simbox)
    acc = newton(phi, pos, simbox)
    vel .+= acc .* Fa(a,Om_l,Om_m) * da/2
    vel ./= simbox.a_f # back to velocities
    if simbox.steps%sp != 0 save(output * "pos_" * string(simbox.steps) * ".jld", "pos", pos, "z", round(simbox.z_f,digits=2)) end
end

"integrate using kick-drift-kick scheme"
function kickNdrift!(pos::Vector{Vector{T}},vel::Vector{Vector{T}}, simbox::SimBox) where{T<:AbstractFloat}

    boxsize = convert(simbox.prec, simbox.boxsize)
    Om_m = simbox.Om_m
    Om_l = simbox.Om_l
    h = simbox.h
    da = (simbox.a_f-simbox.a_i)/simbox.steps
    a  = simbox.a_i
    vel .*= (simbox.a_i-da/2) # rewrite as momenta, eq(6,8) K97, x0H0=1

    for s in 1:simbox.steps
	println("step= ",s)
        phi = poisson(pos, a, simbox)
        acc = newton(phi, pos, simbox)
        vel .+= acc .* Fa(a,Om_l,Om_m) * da        # kick
        a += da/2
	pos .+= vel .* Fa(a,Om_l,Om_m) / a^2 * da  # drift
	a += da/2
	PBC!(pos,boxsize)
    end
    phi = poisson(pos, a, simbox)
    acc = newton(phi, pos, simbox)
    vel .+= acc .* Fa(a,Om_l,Om_m) * da/2
    vel ./= simbox.a_f # back to velocities
end
