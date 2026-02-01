# Code for interpolating the twisted pair automorphism according to Section 8.3.1.1

using DecomposingPolynomialSystems, DecomposingGroupRepresentations
using LinearAlgebra: norm

include("../utils.jl")

# 1. Run monodromy and extract the twisted pair automrophism
include("monodromy.jl")
D = aut_permutations(F)
twisted_pair = D[2]
sol1_image = twisted_pair[1]

# 2. Decompose the vector space of bounded-degree polynomials into isotypic components
include("repr-theory/Rtabxy.jl")
unknwns = Set(vcat(R[:], t, α, β))
μs = useful_highest_weights(vars_irrs, unknwns)
length(μs) # number of useful nullspace computations
max_dim, npairs = tested_hws_pairs_all_μs(iso, μs) # (size of largest Vandermonde matrix, number of all pairs (Hλ₁, Hλ₂))

# 3. Collect samples of the parametric system for the interpolation
sample!(F; path_ids = [1, sol1_image], n_instances = max_dim) # track only 2 solutions that are in correspondece under Ψ
sols = samples(F)[[1, sol1_image]].solutions # dim: nuknowns x nsolutions x ninstances
params = samples(F)[[1, sol1_image]].parameters # dim: nparameters x ninstances

# 4. Interpolate the highest weight vector ψ_R of <Ψ_R> according to Section 8.3.1.1
unknowns(F)[1:9] # rotation matrix unknowns
smpls = vcat(sols[:,1,:], params) # samples corresponding to the solution with idx 1
evals = zeros(ComplexF64, max_dim) # evaluations of ψ_R at smpls
for i in 1:max_dim
    evals[i] = sols[1,2,i] + im*sols[2,2,i] + im*sols[4,2,i] - sols[5,2,i]
end
ν = Weight([1,1,0,0,0,0,0,0,0,0,0,0,0]) # weight of ψ_R
for λ₁ in highest_weights(iso) # for each weight of the numerator of ψ_R
    λ₂ = λ₁ - ν # weight of the denominator of ψ_R
    if has_weight(iso, λ₂)
        Vnum = iso[λ₁] # isotypic component for the numerator
        Vden = iso[λ₂] # isotypic component for the denominator
        A = vandermonde_matrix(Vnum, Vden, vars, smpls, evals) # Vandermonde matrix from (7.7)
        N = Matrix(transpose(nullspace(A; rtol = 1e-10))) # transposed nullspace of A
        DGR.rref!(N, 1e-5) # rref is necessary to identify spurious rows (zero numerator or zero denominator)
        DecomposingGroupRepresentations.sparsify!(N, 1e-5)
        for n in eachrow(N)
            # Interested only in representatives with nonzero numerator or nonzero denominator
            if norm(n[1:mul(Vnum)]) > 1e-8 && norm(n[mul(Vnum)+1:end]) > 1e-8
                println("Nullspace dimension for weights $λ₁ -> $λ₂: ", size(N, 1))
                continue
            end
        end
    end
end # --> The only meaningful weight pair is (λ₁,λ₂) = ([1,1,0,0,0,0,0,0,0,0,0,0,2], [0,0,0,0,0,0,0,0,0,0,0,0,2])

# Focus on the only meaningful weight pair
λ₁ = Weight([1,1,0,0,0,0,0,0,0,0,0,0,2])
λ₂ = λ₁ - ν
has_weight(iso, λ₂)
Vnum, Vden = iso[λ₁], iso[λ₂]
A = vandermonde_matrix(Vnum, Vden, vars, smpls, evals) # Vandermonde matrix from (7.7)
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
