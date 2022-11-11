#######################################
#grid utilities
# HE: Hockney, R.W., & Eastwood, J.W. (1988). Computer Simulation Using Particles (1st ed.)
#######################################
"Account for PBC and avoid floating point error at borders"
function PBC!(pos::Vector{Vector{T}}, boxsize::AbstractFloat) where{T<:AbstractFloat}
    for i in eachindex(pos[1])
          pos[1][i] = mod(pos[1][i],boxsize)
	  if pos[1][i]==boxsize pos[1][i]=0.0 end
          pos[2][i] = mod(pos[2][i],boxsize)
	  if pos[2][i]==boxsize pos[2][i]=0.0 end
          pos[3][i] = mod(pos[3][i],boxsize)
	  if pos[3][i]==boxsize pos[3][i]=0.0 end
    end
end

"Cloud in Cell interpolation"
function CiC(pos::Vector{Vector{T}}, ng::Int, boxsize::AbstractFloat)::Array{T,3} where{T<:AbstractFloat}
    cell = boxsize/ng
    np = size(pos[1])[1]
    posx, posy, posz = pos
    den = zeros(T, (ng,ng,ng))
    for i in 1:np
      cellx = floor(Int, posx[i]/cell)
      celly = floor(Int, posy[i]/cell)
      cellz = floor(Int, posz[i]/cell)
      dx = posx[i]/cell - cellx
      dy = posy[i]/cell - celly
      dz = posz[i]/cell - cellz
      den[1+cellx, 1+celly, 1+cellz] += (1-dx)*(1-dy)*(1-dz)
      den[1+(cellx+1)%ng, 1+celly, 1+cellz] += dx*(1-dy)*(1-dz)
      den[1+cellx, 1+(celly+1)%ng, 1+cellz] += (1-dx)*dy*(1-dz)
      den[1+cellx, 1+celly, 1+(cellz+1)%ng] += (1-dx)*(1-dy)*dz
      den[1+(cellx+1)%ng, 1+(celly+1)%ng, 1+cellz] += dx*dy*(1-dz)
      den[1+(cellx+1)%ng, 1+celly, 1+(cellz+1)%ng] += dx*(1-dy)*dz
      den[1+cellx, 1+(celly+1)%ng, 1+(cellz+1)%ng] += (1-dx)*dy*dz
      den[1+(cellx+1)%ng, 1+(celly+1)%ng, 1+(cellz+1)%ng] += dx*dy*dz
    end

    # transform it to overdensity
    mden = 0.0
    for i in eachindex(den)
      mden += den[i]
    end
    mden /= ng^3

    for i in eachindex(den)
      den[i] /= mden
      den[i] -= 1.0
    end
    return den
end

"assign from mesh to particles using trilinear interpolation"
function trilinear(den::Array{T, 3}, pos::Vector{Vector{T}}, boxsize::AbstractFloat)::Array{T} where{T<:AbstractFloat}
    np = size(pos[1])[1] # number of particles
    posx, posy, posz = pos
    ng = size(den)[1]    # size of mesh 
    cell = boxsize/ng
    mass = zeros(T, np)
    for i in 1:np
        cellx = floor(Int, posx[i]/cell)
        celly = floor(Int, posy[i]/cell)
        cellz = floor(Int, posz[i]/cell)
        dx = posx[i]/cell - cellx
        dy = posy[i]/cell - celly
        dz = posz[i]/cell - cellz
        mass[i] += den[1+cellx, 1+celly, 1+cellz] * (1-dx)*(1-dy)*(1-dz)
        mass[i] += den[1+(cellx+1)%ng, 1+celly, 1+cellz] * dx*(1-dy)*(1-dz)
        mass[i] += den[1+cellx, 1+(celly+1)%ng, 1+cellz] * (1-dx)*dy*(1-dz)
        mass[i] += den[1+cellx, 1+celly, 1+(cellz+1)%ng] * (1-dx)*(1-dy)*dz
        mass[i] += den[1+(cellx+1)%ng, 1+(celly+1)%ng, 1+cellz] * dx*dy*(1-dz)
        mass[i] += den[1+(cellx+1)%ng, 1+celly, 1+(cellz+1)%ng] * dx*(1-dy)*dz
        mass[i] += den[1+cellx, 1+(celly+1)%ng, 1+(cellz+1)%ng] * (1-dx)*dy*dz
        mass[i] += den[1+(cellx+1)%ng, 1+(celly+1)%ng, 1+(cellz+1)%ng] * dx*dy*dz
    end
    return mass
end
