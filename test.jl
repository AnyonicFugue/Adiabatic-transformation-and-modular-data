include("utils.jl")
include("adiabatic.jl")

using LinearAlgebra
import Combinatorics
import SkewLinearAlgebra
import PyPlot


n_step=32
n_length=3
debug=false


function BerryPhase_square(nx::Int,ny::Int,xPBC::Bool,yPBC::Bool)
   A0=construct_A_square(nx,ny,-2.0,1.0,1.0,xPBC,yPBC,0,0) # No x-shift
   A1=construct_A_square(nx,ny,-2.0,1.0,1.0,xPBC,yPBC,1,0) # x-shift for length 1
   res,phase_arr,overlap_with_A0=adiabatic_ver1(A0,A1,n_step,debug)
   display(res)
   display(overlap_with_A0)
   display(phase_arr)
   return res
   # display(phase_arr)
   # display(overlap_with_A0)
end

# Use Pyplot to plot the interactions in the honeycomb lattice.
# The input is the A matrix.

phase_arr=zeros(Float64,n_length)
N_arr=zeros(Int64,n_length)

for l in range(3,3+n_length-1)
   phase_arr[l-2]=BerryPhase_square(l,l,true,true)
   N_arr[l-2]=l^2
end

# Perform linear fit to phase_arr with respect to l^2.
Fit_and_plot(phase_arr,N_arr)

