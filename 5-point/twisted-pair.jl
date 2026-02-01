# Code for interpolating the twisted pair automorphism according to Section 8.3.1.1

include("repr-theory/Rtabxy.jl")
include("monodromy.jl")

max_dim # size of largest Vandermonde matrix
sample!(F; path_ids = [1, sol1_image], n_instances = max_dim)

sols = samples(F)[[1, sol1_image]].solutions
params = samples(F)[[1, sol1_image]].parameters

# Interpolate the highest weight vector of <Ψ_R>
unknowns(F)[1:9]
evals = zeros(ComplexF64, max_dim)
for i in 1:max_dim
    evals[i] = sols[1,2,i] + im*sols[2,2,i] + im*sols[4,2,i] - sols[5,2,i]
end
smpls = vcat(sols[:,1,:], params)
ν = Weight([1,1,0,0,0,0,0,0,0,0,0,0,0])
for λ₁ in highest_weights(iso)
    λ₂ = λ₁ - ν
    if has_weight(iso, λ₂)
        V₁ = iso[λ₁]
        V₂ = iso[λ₂]
        A = vandermonde_matrix(V₁, V₂, vars, smpls, evals)
        println(size(A))
        N = nullspace(A; rtol = 1e-10)
        DecomposingGroupRepresentations.sparsify!(N, 1e-5)
        println("Nullspace dimension for weights $λ₁ -> $λ₂: ", size(N, 2))
    end
end

λ₁ = Weight([1,1,0,0,0,0,0,0,0,0,0,0,2])
λ₂ = λ₁ - ν
has_weight(iso, λ₂)
Vnum, Vden = iso[λ₁], iso[λ₂]
A = vandermonde_matrix(Vnum, Vden, vars, smpls, evals)
N = nullspace(A)

hw_vectors(Vnum; as_vectors = true)
hw_vectors(Vden; as_vectors = true)

# ----------------------------------------------------------------------
λ₁ = Weight([1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 2])
λ₂ = λ₁ - ν
has_weight(iso, λ₂)
Vnum, Vden = iso[λ₁], iso[λ₂]
A = vandermonde_matrix(Vnum, Vden, vars, smpls, evals)
N = nullspace(A)

hw_vectors(Vnum; as_vectors = true)
hw_vectors(Vden; as_vectors = true)