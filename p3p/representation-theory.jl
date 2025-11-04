using DecomposingGroupRepresentations

@polyvar R[1:3,1:3] t[1:3] α[1:3] x[1:3,1:3] X[1:4,1:3] n[1:3]

SO3 = LieGroup("SO", 3)
a₁ = MatrixGroupAction(SO3, vcat(eachcol(R), [t], eachcol(x)))
a₂ = MatrixGroupAction(SO3, vcat(eachrow(R), eachcol(X[1:3,:])))

E = zeros(Int, 7, 30) # α x X t n
[(E[i,i], E[i,(1+3*i):(3+3*i)]) = (1,[-1,-1,-1]) for i in 1:3]
[(E[3+i,(1+3*i):(3+3*i)], E[3+i, (9+4*i):(12+4*i)]) = (ones(3), ones(4)) for i in 1:3]
(E[7,25:27], E[7,[16,20,24,28,29,30]]) = (ones(3), -ones(6))
T = ScalingLieGroup(E)
a₃ = ScalingLieGroupAction(T, vcat(α, x[:], X[:], t, n))

A = a₁ × a₂ × a₃

vars = vcat(R[:], t, α, x[:], X[:], n)
V = BoundedDegreePolynomials(vars, 3)
ρ = GroupRepresentation(A, V)
irrs = irreducibles(ρ)
iso = isotypics(irrs)
