using LinearAlgebra

# hbar=1, m=1

################### Potential #######################

V(x::Float64) = 1000x.^2



################### Domain ##########################

N = 1e3
a,b = 0,1
x = Vector(a:(b-a)/N:b)
h = x[2]-x[1]

################################### Matrices ###################################

function crear_T(N)
    M = zeros(Float64, N-2, N-2)
    for i in 1:1:N-2
        for j in 1:1:N-2
            if abs(i-j)==1
                M[j,i]=1
            elseif i==j
                M[j,i]=-2
            end
        end
    end
    M
end

T(h) = (1/(2(h)^2))*crear_T(int)


function crear_V(N)
    M = zeros(Float64, N-2, N-2)
    for i in 1:1:N-2
        for j in 1:1:N-2
            if i==j
                M[j,i]=V(x[i])
            end
        end
    end
    M
end



################ Hamiltonian ###########################

ℍ = (1/(2(h)^2))*crear_T(Int(N)) + crear_V(Int(N))

auto_vect = eigvecs(ℍ)

auto_vals = eigvals(ℍ)

Energias =	sort!(auto_vals)

Energias = Energias/minimum(auto_vals)
################ Plot ################################

using Plots
n_plots=10
plt = plot()
for i in 1:4
    eigen_vect = [0.]
    eigen_vect=append!(append!(eigen_vect,auto_vect[:,i]),[0.])
    plot!(eigen_vect,label="ψ$i",markersize = 1)
    display(plt)
end
################## Amplitude ###############################
using Plots
n_plots=10
plt = plot()
for i in 1:4
    eigen_vect = [0.]
    eigen_vect=append!(append!(eigen_vect,auto_vect[:,i]),[0.])
    plot!(eigen_vect.^2,label="|ψ$i|²",markersize = 1)
    display(plt)
end
