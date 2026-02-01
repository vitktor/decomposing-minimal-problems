# Code for interpolating the twisted pair automorphism according to Section 8.3.1.1

using DecomposingPolynomialSystems, DecomposingGroupRepresentations
using LinearAlgebra: norm

include("../utils.jl")
include("monodromy.jl")
include("repr-theory/Rtabxy.jl")

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
        Vnum = iso[λ₁]
        Vden = iso[λ₂]
        A = vandermonde_matrix(Vnum, Vden, vars, smpls, evals)
        N = Matrix(transpose(nullspace(A; rtol = 1e-10)))
        DGR.rref!(N, 1e-5)
        DecomposingGroupRepresentations.sparsify!(N, 1e-5)
        for n in eachrow(N)
            # Interested only in representatives with nonzero numerator or nonzero denominator
            if norm(n[1:mul(Vnum)]) > 1e-8 && norm(n[mul(Vnum)+1:end]) > 1e-8
                println("Nullspace dimension for weights $λ₁ -> $λ₂: ", size(N, 2))
                continue
            end
        end
    end
end # --> The only meaningful weight pair is (λ₁,λ₂) = ([1,1,0,0,0,0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,0,0,0,0,0,2])

# Focus on the only meaningful weight λ₁ = [1,1,0,0,0,0,0,0,0,0,0,0,2]
λ₁ = Weight([1,1,0,0,0,0,0,0,0,0,0,0,2])
λ₂ = λ₁ - ν
has_weight(iso, λ₂)
Vnum, Vden = iso[λ₁], iso[λ₂]
A = vandermonde_matrix(Vnum, Vden, vars, smpls, evals)
N = nullspace(A)
n = N[:,1]
n = n / n[argmax(abs.(n))]
DecomposingGroupRepresentations.sparsify!(n, 1e-10)

# Reconstruct numerator and denominator of ψ_R from the nullspace coefficients
hwv_num = hw_vectors(Vnum; as_vectors = true)
hwv_den = hw_vectors(Vden; as_vectors = true)
a = div_by_smallest_coeff(sum(n[1:mul(Vnum)].*hwv_num))
b = div_by_smallest_coeff(sum(n[mul(Vnum)+1:end].*hwv_den))

# Verify that ψ_R = a/b corresponds to the well-known formula for the twisted pair
M = (2*t*t' - t'*t*I(3)) * R
f = M[1,1] + im*M[1,2] + im*M[2,1] - M[2,2]
g = t'*t
f/g == a/b
f == -a || f == a
g == -b || g == b
