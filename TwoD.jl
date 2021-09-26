#######################################################################

# Naive 2D NN solver___Version 0.1
# by: Jose Miguel Muñoz A
########################################################################

########################################################################3
#                   Packages
########################################################################
using LinearAlgebra, BenchmarkTools, Flux, Plots,Interpolations, Random
using Distributions
import PlotlyJS
using SparseArrays
using Kronecker
using Arpack

#    ħ=1, m=1
########################################################################
#                   Space discretization
########################################################################

# Mesh size
Nx =  110
Ny = Nx
# Range information
ax,bx = -1,1
ay,by = -1,1
# Space grids
Rx = ax:(bx-ax)/(Nx):bx
Ry = ay:(by-ay)/(Nx):by

########################################################################
#                   POTENTIAL
########################################################################
# Define the potential function
V(x,y) =  (x^2+y^2)

"""
Map the potential to each point in space
"""
function make_vs()
   z = zeros(Nx,Ny)
   for px in 1:Nx
      for py in 1:Ny
         z[px,py]=V(Rx[px],Ry[py])
   end end
   z
end

𝓥 = make_vs()
# Plot the potential
PlotlyJS.plot(PlotlyJS.contour(x=Rx, y=Ry,z=𝓥, contours_coloring="heatmap" ))


########################################################################
#              Construct the strange matrices
########################################################################

#Create the diagonal sparse matrix to use
D = spdiagm(-1 => ones(Nx-1), 0 => ones(Nx)*(-2) ,1 => ones(Nx-1) )

# Define the kinetic energy as the Kronecker sum
𝕋 = sparse((D ⊕ D)).*(-0.5)
# Define the potential function
𝕌 = spdiagm(0 => reshape((𝓥),Nx^2))
# Create the Hamiltonian
ℍ = 𝕋 + 𝕌
########################################################################
#              Get the eigensystem
########################################################################
# Use an approximation for the 20 first eigensystems of the problem
λ, ϕ = eigs(ℍ, nev = 20, which=:SM)
# Mat the resulting function to the initial space
Φ(n) =  reshape((ϕ)[:,n],Nx,Nx)

########################################################################
#              See the ψ
########################################################################
PlotlyJS.plot(PlotlyJS.heatmap(x=Rx, y=Ry,z=Φ(1), showscale=false, connectgaps=true, zsmooth="best" ))

########################################################################
#              See the ψ²
########################################################################
PlotlyJS.plot(PlotlyJS.heatmap(x=Rx, y=Ry,z=Φ(4).^2, showscale=false, connectgaps=true, zsmooth="best" ))

########################################################################
#              Study the linear relation
########################################################################
# Evaluate the energy levels by normalization
ϐ = λ[1]/2
ε = λ/ϐ
# Plot the energy levels
Plots.scatter(ε)

########################################################################
#             Some Flexing
########################################################################
# Create the old fasion animation
gr()
using Plots
anim = @animate for n = 1:15
    heatmap(Φ(n).^2)
end
gif(anim, fps = 3)
