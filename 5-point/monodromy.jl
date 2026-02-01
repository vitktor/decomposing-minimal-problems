using DecomposingPolynomialSystems
using LinearAlgebra: nullspace, I

@var R[1:3,1:3] t[1:3] α[1:5] β[1:5] x[1:3,1:5] y[1:3,1:5] a[1:4]
eqs = vcat((R'*R - I)[:], [det(R) - 1])
for i in 1:5
    append!(eqs, β[i]*y[:,i] - R*α[i]*x[:,i] - t)
end
push!(eqs, a'*[t; 1])

F = ParametricSystem(eqs; unknowns = vcat(R[:], t, α, β), parameters = vcat(x[:], y[:], a))

function fabricateSample()
    R₁ = Matrix{ComplexF64}(I(3))
    t₁ = zeros(ComplexF64, 3)
    R₂ = c2R(randn(ComplexF64, 3))
    t₂ = randn(ComplexF64, 3)
    X = randn(ComplexF64, 3, 5)
    x = [R₁ t₁]*a2p(X)
    y = [R₂ t₂]*a2p(X)
    α, β = ones(ComplexF64, 5), ones(ComplexF64, 5)
    n = nullspace(reshape([t₂; 1], 1, 4))
    a = vec(n*randn(ComplexF64, size(n, 2)))
    return (vcat(R₂[:], t₂, α, β), vcat(x[:], y[:], a))
end

x₀, p₀ = fabricateSample()

F = run_monodromy(F, ([x₀], p₀))
