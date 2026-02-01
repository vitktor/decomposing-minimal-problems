import DynamicPolynomials as DP
import HomotopyContinuation as HC
to_expression(
    f::DP.AbstractPolynomialLike,
    vars_DP::Vector{<:DP.Variable},
    vars_HC::Vector{<:HC.Variable}
) = Expression(DP.subs(f, vars_DP => vars_HC))

f = hw_vectors(iso[first(highest_weights(iso))]; as_vectors = true)[1]
vars_HC = vcat(unknowns(F), parameters(F))
vars
to_expression(f, vars, vars_HC)

