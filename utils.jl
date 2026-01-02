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
    unknowns::Set{<:Variable}
)
    μs = Weight[]
    for irr in vars_irrs
        hwv = hw_vector(irr)
        if variables(vector(hwv)) ⊆ unknowns
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