include("utils.jl")
using LinearAlgebra


A=construct_A(10,10,1.0,1.1,1.2,true,true)
# println(Hermitian(1im*A))
# println(A)



e2=eigen(Hermitian(1im*A))

println("values")
println(e2.values)
println("vectors")
# println(e2.vectors)