include("utils.jl")
using LinearAlgebra
import Combinatorics
import SkewLinearAlgebra
import PyPlot

A=construct_A_square(3,3,-2.0,1.0,1.0,true,true,0,0)
display(A[1:8,1:8])

# Use Pyplot to plot the interactions in the honeycomb lattice.
# The input is the A matrix.


# plot_A(A)