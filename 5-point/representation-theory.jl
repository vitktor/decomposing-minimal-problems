using DecomposingGroupRepresentations

@polyvar R[1:3,1:3] t[1:3] α[1:5] β[1:5] x[1:3,1:5] y[1:3,1:5]

SO3 = LieGroup("SO", 3)
a₁ = MatrixGroupAction(SO3, vcat(eachcol(R), [t], eachcol(y)))
a₂ = MatrixGroupAction(SO3, vcat(eachrow(R), eachcol(x)))

E = zeros(Int, 11, 43) # α β x y t
E[11, vcat(1:10, 41:43)] .= 1
[E[i,i] = -1 for i in 1:10]
[E[i, (11+3*(i-1)):(13+3*(i-1))] .= 1 for i in 1:10]
T = ScalingLieGroup(E)
a₃ = ScalingLieGroupAction(T, vcat(α, β, x[:], y[:], t))

A = a₁ × a₂ × a₃

vars = vcat(R[:], t, α, β, x[:], y[:])
V = BoundedDegreePolynomials(vars, 3)
ρ = GroupRepresentation(A, V)
irrs = irreducibles(ρ)
iso = isotypics(irrs)