using DecomposingGroupRepresentations

@polyvar E[1:3,1:3] x[1:3,1:5] y[1:3,1:5]

SO3 = LieGroup("SO", 3)
a₁ = MatrixGroupAction(SO3, vcat(eachcol(E), eachcol(y)))
a₂ = MatrixGroupAction(SO3, vcat(eachrow(E), eachcol(x)))

M = zeros(Int, 11, 39) # E x y
M[11, 1:9] .= 1
[M[i, (10+3*(i-1)):(12+3*(i-1))] .= 1 for i in 1:10]
T = ScalingLieGroup(M)
a₃ = ScalingLieGroupAction(T, vcat(E[:], x[:], y[:]))

A₁ = a₁ × a₂ × a₃
A₂ = a₃
A = A₁ # choice of a group action

vars = vcat(E[:], x[:], y[:])
V = VariableSpace(vars)
ρ = GroupRepresentation(A, V)
vars_irrs = irreducibles(ρ)

D = 3
V = BoundedDegreePolynomials(vars, D)
ρ = GroupRepresentation(A, V)
irrs = irreducibles(ρ)
iso = isotypics(irrs)

iso[Weight([1,1,0,0,0,0,0,0,0,0,0,0,3])]
iso[Weight([0,0,0,0,0,0,0,0,0,0,0,0,3])]
basis(iso[Weight([0,0,1,0,0,0,0,1,0,0,0,0,1])])