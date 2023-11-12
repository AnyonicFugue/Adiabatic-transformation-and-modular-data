include("utils.jl")
using LinearAlgebra
import Combinatorics
import SkewLinearAlgebra

m=rand(Float64,6,6)
SkewLinearAlgebra.skewhermitian!(m)

println(Pfaffian(m))
println(SkewLinearAlgebra.pfaffian(m))


