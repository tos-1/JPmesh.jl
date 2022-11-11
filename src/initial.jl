###################################################
#generate initial conditions with 2LPT
#B01: https://arxiv.org/abs/astro-ph/0112551
#EHU97: https://arxiv.org/abs/astro-ph/9710252
###################################################

"It returns the density field on the Fourier grid, with power Pk=Ak^-ns"
function Gaussimple(simbox::SimBox)::Array{Complex{SimBox().prec}, 3}
    ng = simbox.ng
    boxsize = simbox.boxsize
    k = kgrid(ng,boxsize)
    pk = Array{simbox.prec,3}(undef, div(ng,2)+1, ng, ng)

    kmin = 2 * pi / boxsize
    for i in 1:div(ng,2)+1, j in 1:ng, k in 1:ng
          kx = if i<=div(ng,2)+1 i-1 else (i-1-ng) end
          ky = if j<=div(ng,2)+1 j-1 else (j-1-ng) end
          kz = if k<=div(ng,2)+1 k-1 else (k-1-ng) end
	  ksum = sqrt(kx^2 + ky^2 + kz^2)*kmin
	  pk[i,j,k] = simbox.s8 * (ksum / kmin)^(-simbox.ns)
    end
    pk[1,1,1] = 0.0

    # generate Gaussian white noise
    rng = MersenneTwister(simbox.seed)
    d = randn(rng, simbox.prec, (ng, ng, ng))
    dk = rfft(d)
    dk .*= sqrt.(pk) .* ng^1.5 / boxsize^1.5
    rfft(irfft(dk,ng)) # make hermitian
end

"It returns the density field on the Fourier grid, Pk from EHU97.
 Already correctly normalized at simbox.z with simbox.s8"
function Gauss(simbox::SimBox)::Array{Complex{SimBox().prec}, 3}
    ng = simbox.ng
    boxsize = simbox.boxsize
    pk = EisHu(simbox)

    # generate Gaussian white noise
    rng = MersenneTwister(simbox.seed)
    d = randn(rng, simbox.prec, (ng, ng, ng)) # define precision here!
    dk = rfft(d)

    kmin = 2 * pi / boxsize
    for i in 1:div(ng,2)+1, j in 1:ng, k in 1:ng
          kx = if i<=div(ng,2)+1 i-1 else (i-1-ng) end
          ky = if j<=div(ng,2)+1 j-1 else (j-1-ng) end
          kz = if k<=div(ng,2)+1 k-1 else (k-1-ng) end
	  ksum = sqrt(kx^2 + ky^2 + kz^2)*kmin
	  dk[i,j,k] *= sqrt(pk(ksum)) * ng^1.5 / boxsize^1.5
    end
    dk[1,1,1] = 0.0
    rfft(irfft(dk,ng)) # make hermitian
end

"Anti-divergence operator in Fourier space. Divergence in Fourier space as input."
function invdiv(psik::Array{T, 3}, simbox::SimBox) where{T<:Complex}
    ng = simbox.ng
    boxsize = simbox.boxsize

    kmin = 2 * pi / boxsize
    phicx = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    phicy = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    phicz = Array{T, 3}(undef, div(ng,2)+1, ng, ng)

    for i in 1:div(ng,2)+1, j in 1:ng, k in 1:ng
        kx = if i<=div(ng,2)+1 i-1 else (i-1-ng) end
        ky = if j<=div(ng,2)+1 j-1 else (j-1-ng) end
        kz = if k<=div(ng,2)+1 k-1 else (k-1-ng) end
        kk = sqrt(kx^2+ky^2+kz^2) * kmin
        phicx[i,j,k] = -1im * kx * kmin * psik[i,j,k]/kk^2
        phicy[i,j,k] = -1im * ky * kmin * psik[i,j,k]/kk^2
        phicz[i,j,k] = -1im * kz * kmin * psik[i,j,k]/kk^2
    end
    # singularity in kk=0
    phicx[1,1,1] = 0.0
    phicy[1,1,1] = 0.0
    phicz[1,1,1] = 0.0

    psix = irfft(phicx, ng)
    psiy = irfft(phicy, ng)
    psiz = irfft(phicz, ng)

    return [psix, psiy, psiz]
end

"second order displacemet potential"
function twolpt(dk::Array{T, 3},simbox::SimBox)::Array{T,3} where{T<:Complex}
    ng = simbox.ng
    boxsize = simbox.boxsize

    phicxx = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    phicxy = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    phicxz = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    phicyy = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    phicyz = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    phiczz = Array{T, 3}(undef, div(ng,2)+1, ng, ng)
    for i in 1:div(ng,2)+1, j in 1:ng, k in 1:ng
        kx = if i<=div(ng,2)+1 i-1 else (i-1-ng) end
        ky = if j<=div(ng,2)+1 j-1 else (j-1-ng) end
        kz = if k<=div(ng,2)+1 k-1 else (k-1-ng) end
        kk = sqrt(kx^2+ky^2+kz^2)
        phicxx[i,j,k] = kx * kx * dk[i,j,k]/kk^2
        phicxy[i,j,k] = kx * ky * dk[i,j,k]/kk^2
        phicxz[i,j,k] = kx * kz * dk[i,j,k]/kk^2
        phicyy[i,j,k] = ky * ky * dk[i,j,k]/kk^2
        phicyz[i,j,k] = ky * kz * dk[i,j,k]/kk^2
        phiczz[i,j,k] = kz * kz * dk[i,j,k]/kk^2
    end
    # singularity in kk=0
    phicxx[1,1,1] = 0.0
    phicxy[1,1,1] = 0.0
    phicxz[1,1,1] = 0.0
    phicyy[1,1,1] = 0.0
    phicyz[1,1,1] = 0.0
    phiczz[1,1,1] = 0.0

    phirxx = irfft(phicxx, ng)
    phirxy = irfft(phicxy, ng)
    phirxz = irfft(phicxz, ng)
    phiryy = irfft(phicyy, ng)
    phiryz = irfft(phicyz, ng)
    phirzz = irfft(phiczz, ng)

    # eq.(103) B01
    phi2 = (phirxx .* phiryy .+ phirxx .* phirzz .+ phiryy .* phirzz .- 
		phirxy .* phirxy .- phirxz .* phirxz .- phiryz .* phiryz)

    rfft(phi2)
end


"It returns the initial conditions via 2LPT"
function InitialConditions(simbox::SimBox)::Tuple{ Vector{Vector{SimBox().prec}} , Vector{Vector{SimBox().prec}} }
    Om_m = simbox.Om_m 
    Om_l = simbox.Om_l

    a_i = simbox.a_i
    da = (simbox.a_f-simbox.a_i)/simbox.steps
    a_m1 = a_i - da/2.0
    @assert a_m1 > 0

    ng = simbox.ng
    boxsize = simbox.boxsize
    cell = boxsize/ng

    # zeldovich approx
    dk = Gauss(simbox)

    # eq.(97) B01
    D2(a) =  - 3. / 7. * Om0z(a_i, Om_l, Om_m)^(-1/143.) * D1(a, Om_m, Om_l)^2.
    psik_2lpt = twolpt(dk, simbox )

    psi      = invdiv(-dk, simbox)       # za at a_i
    psi_2lpt = invdiv(psik_2lpt, simbox)
    psik_za  = Nothing
    psik_2lpt = Nothing
    GC.gc()

    # eq.(101) B01
    f1 = Om0z(a_m1, Om_l, Om_m)^(5/9)
    f2 = 2 * Om0z(a_m1, Om_l, Om_m)^(6/11)

    # eq.(99) B01, in units of Mpc/(h *H0 * time), it starts at n-1/2
    vel = f1 * (D1(a_m1, Om_m, Om_l)/D1(a_i, Om_m, Om_l)) * a_m1 * E(a_m1,Om_l,Om_m) * psi + f2 * (D2(a_m1)/D1(a_i, Om_m, Om_l)^2) * a_m1 * E(a_m1,Om_l,Om_m) * psi_2lpt

    # total displacement
    psi .+= psi_2lpt * (D2(a_i)/D1(a_i, Om_m, Om_l)^2.)
    psi_2lpt = Nothing
    GC.gc()

    # initial positions
    for i in 1:ng
      for j in 1:ng
        for k in 1:ng
          psi[1][i,j,k] = psi[1][i,j,k] + (i-1)*cell
          psi[2][i,j,k] = psi[2][i,j,k] + (j-1)*cell
          psi[3][i,j,k] = psi[3][i,j,k] + (k-1)*cell
        end
      end
    end

    pos = [ reshape(psi[1],ng^3), reshape(psi[2],ng^3), reshape(psi[3],ng^3) ]
    psi = Nothing
    GC.gc()
    PBC!(pos, boxsize)
    vel = [ reshape(vel[1],ng^3), reshape(vel[2],ng^3), reshape(vel[3],ng^3) ]
    return pos, vel
end

#simbox = SimBox(ng=128, boxsize=250.0, z=0)
#dk=Gauss(simbox)
#den=irfft(dk, 128)
#display( heatmap( 1:size(den,1), 1:size(den,2), den[:,:,2]) )
#gui()
#pos, vel = InitialConditions(simbox)
#den = CiC(pos, simbox.ng, simbox.boxsize);
#display( heatmap( 1:size(den,1), 1:size(den,2), den[:,:,2]) )
#gui()
#display( heatmap( 1:size(d,1), 1:size(d,2), d[:,:,32]) )
#psi,vel = InitialConditions(0.001,0.7,0.3,32,100.,10.0,-3.0,12)
#println(typeof(psi), axes(psi))
#println(typeof(psi[1]), axes(psi[1]))

