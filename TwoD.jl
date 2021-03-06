#######################################################################

# Naive 2D NN solver___Version 0.1
# by: Jose Miguel MuÃ±oz A
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

#    Ä§=1, m=1
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

ð¥ = make_vs()
# Plot the potential
PlotlyJS.plot(PlotlyJS.contour(x=Rx, y=Ry,z=ð¥, contours_coloring="heatmap" ))


########################################################################
#              Construct the strange matrices
########################################################################

#Create the diagonal sparse matrix to use
D = spdiagm(-1 => ones(Nx-1), 0 => ones(Nx)*(-2) ,1 => ones(Nx-1) )

# Define the kinetic energy as the Kronecker sum
ð = sparse((D â D)).*(-0.5)
# Define the potential function
ð = spdiagm(0 => reshape((ð¥),Nx^2))
# Create the Hamiltonian
â = ð + ð
########################################################################
#              Get the eigensystem
########################################################################
# Use an approximation for the 20 first eigensystems of the problem
Î», Ï = eigs(â, nev = 20, which=:SM)
# Mat the resulting function to the initial space
Î¦(n) =  reshape((Ï)[:,n],Nx,Nx)

########################################################################
#              See the Ï
########################################################################
PlotlyJS.plot(PlotlyJS.heatmap(x=Rx, y=Ry,z=Î¦(1), showscale=false, connectgaps=true, zsmooth="best" ))

########################################################################
#              See the ÏÂ²
########################################################################
PlotlyJS.plot(PlotlyJS.heatmap(x=Rx, y=Ry,z=Î¦(4).^2, showscale=false, connectgaps=true, zsmooth="best" ))

########################################################################
#              Study the linear relation
########################################################################
# Evaluate the energy levels by normalization
Ï = Î»[1]/2
Îµ = Î»/Ï
# Plot the energy levels
Plots.scatter(Îµ)

########################################################################
#             Some Flexing
########################################################################
# Create the old fasion animation
gr()
using Plots
anim = @animate for n = 1:15
    heatmap(Î¦(n).^2)
end
gif(anim, fps = 3)
