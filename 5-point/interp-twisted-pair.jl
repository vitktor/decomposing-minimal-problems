include("representation-theory.jl")

M = (2*t*t'-t'*t*eye(3))*R
M[1,1]+im*M[2,1]+im*M[1,2]-M[2,2]