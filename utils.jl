import DynamicPolynomials as DP
import HomotopyContinuation as HC
import DecomposingGroupRepresentations as DGR

"""
Quaternion to scaled rotation matrix.
"""
function q2sR(q)
    a,b,c,d = q
    return [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
            2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
            2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2]
end

"""
Quaternion to rotation matrix.
"""
function q2R(q)
    a,b,c,d = q
    return 1/(a^2+b^2+c^2+d^2)*q2sR(q)
end

"""
Cayley vector to scaled rotation matrix.
"""
function c2sR(c)
    x, y, z = c
    return [1+x^2-y^2-z^2 2*(x*y-z) 2*(y+x*z);
            2*(x*y+z) 1-x^2+y^2-z^2 2*(y*z-x);
            2*(x*z-y) 2*(x+y*z) 1-x^2-y^2+z^2]
end

"""
Cayley vector to rotation matrix.
"""
function c2R(c)
    x, y, z = c
    return 1/(1+x^2+y^2+z^2)*c2sR(c)
end

"""
Quaternion to Cayley vector.
"""
q2c(q) = [q[2]/q[1], q[3]/q[1], q[4]/q[1]]

"""
Rotation matrix to Cayley vector.
"""
R2c(R) = 1/(1+tr(R))*xx2v((R-transpose(R)))

"""
Returns 
"""
function tested_hws_pairs(
    iso::IsotypicDecomposition,
    μ::Weight
)
    max_dim, npairs = 0, 0
    for ic in iso
        λ₂ = highest_weight(ic)
        λ₁ = λ₂ + μ
        if has_weight(iso, λ₁)
            npairs += 1
            new_dim = mul(iso[λ₁])+mul(ic)
            if new_dim > max_dim
                max_dim = new_dim
            end
        end
    end
    return max_dim, npairs
end

function useful_highest_weights(
    vars_irrs::Vector{<:IrreducibleRepresentation},
    unknowns::Set{<:DP.Variable}
)
    μs = Weight[]
    for irr in vars_irrs
        hwv = hw_vector(irr)
        if DGR.variables(vector(hwv)) ⊆ unknowns
            push!(μs, highest_weight(irr))
        end
    end
    return μs
end

function tested_hws_pairs_all_μs(
    iso::IsotypicDecomposition,
    μs::Vector{<:Weight}
)
    max_dim, npairs = 0, 0
    for μ in μs
        d, n = tested_hws_pairs(iso, μ)
        if d > max_dim
            max_dim = d
        end
        npairs += n
    end
    return max_dim, npairs
end

function vandermonde_matrix(
    V_num::IsotypicComponent,
    V_den::IsotypicComponent,
    vars_DP::Vector{<:DP.Variable},
    samples::Matrix{ComplexF64},
    evals::Vector{ComplexF64}
)
    hwv₁ = hw_vectors(V_num; as_vectors = true)
    hwv₂ = hw_vectors(V_den; as_vectors = true)
    k = length(hwv₁) + length(hwv₂)
    A = zeros(ComplexF64, k, k)
    for i in 1:k # for each row of A
        for (j, hwv) in enumerate(hwv₁)
            A[i, j] = hwv(vars_DP => samples[:, i])
        end
        for (j, hwv) in enumerate(hwv₂)
            A[i, j + length(hwv₁)] = -evals[i]*hwv(vars_DP => samples[:, i])
        end
    end
    return A
end

to_expression(
    f,
    vars_DP::Vector{<:DP.Variable},
    vars_HC::Vector{<:HC.Variable}
) = HC.Expression(DP.subs(f, vars_DP => vars_HC))