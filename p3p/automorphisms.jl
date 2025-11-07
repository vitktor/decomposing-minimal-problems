using DecomposingPolynomialSystems

@var R[1:3,1:3] t[1:3] α[1:3] x[1:3,1:3] X[1:4,1:3] n[1:3]
eqs = (R'R-eye(Int, 3))[:] ∪ [det(R)-1] ∪ []