# create gif of evolution

push!(LOAD_PATH, "./src/")
include("src/grid.jl")
using LaTeXStrings
using JLD
using Plots

ng = 128
boxsize = 256.0

# reorder
zs = zeros(0)
order = zeros(0)
files = readdir("data/")
for fil in readdir("data/")
  a = lstrip(fil,['p','o','s','_'])
  b = rstrip(a,['.','j','l','d'])
  append!( order, parse(Int, b))
end
order = sortperm(order)
ordered_files = String[]
for i in order
  push!(ordered_files,files[i])
end

anim = @animate for fil in ordered_files
  pos = load("data/"*fil,"pos")
  z = load("data/"*fil,"z")
  den = CiC(pos, ng, boxsize)
  den = log.(2.0 .+ den)
  heatmap(range(0.,boxsize-boxsize/ng,ng),
	  range(0.,boxsize-boxsize/ng,ng),
	  den[:,:,div(ng,2)], 
	  xlabel=L"Mpc/h",
	  ylabel=L"Mpc/h",
	  title=L"z="*"$(z)",
          colorbar_title = L"\log(2+\delta)",
)
  plot!(size=(600,600))
end

gif(anim, "images/LSS.gif", fps = 1)

