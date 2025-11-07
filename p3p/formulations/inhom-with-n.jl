using DecomposingPolynomialSystems

@var R[1:3,1:3] t[1:3] α[1:3] n[1:3] x[1:3,1:3] X[1:3,1:3]

eqs = vcat(reshape(R'*R - eye(Int,3), 9), [det(R) - 1])
for i in 1:3
    append!(eqs, α[i]*x[:,i] - [R t]*[X[:,i]; 1])
end
append!(eqs, (n'*X)[:] - ones(Int, 3))
F = System(eqs; variables = vcat(R[:], t, α, n), parameters = vcat(X[:], x[:]))

scaling_symmetries(F)


using DecomposingGroupRepresentations

@polyvar R[1:3,1:3] t[1:3] α[1:3] n[1:3] x[1:3,1:3] X[1:3,1:3]

SO3 = LieGroup("SO", 3)
a₁ = MatrixGroupAction(SO3, vcat(eachcol(R), [t], eachcol(x)))
a₂ = MatrixGroupAction(SO3, vcat(eachrow(R), eachcol(X)))

E = zeros(Int, 4, 27) # α x X t n
[E[i,i] = 1 for i in 1:3]
[E[i, (1+3*i):(3+3*i)] .= -1 for i in 1:3]
E[4, vcat(1:3, 13:24)] .= 1
E[4, 25:27] .= -1
T = ScalingLieGroup(E)
a₃ = ScalingLieGroupAction(T, vcat(α, x[:], X[:], t, n))

A = a₁ × a₂ × a₃

vars = vcat(R[:], t, α, x[:], X[:], n)
V = BoundedDegreePolynomials(vars, 3)
ρ = GroupRepresentation(A, V)
irrs = irreducibles(ρ)
iso = isotypics(irrs) # (9139, 450, 32, 288)
