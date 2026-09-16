## Description #############################################################################
#
# Mathematical helpers shared by the celestial body algorithms.
#
############################################################################################

"""
    _evalpoly_with_derivative(x::Number, coefficients::NTuple{N, Number}) -> Number, Number

Evaluate the polynomial with the `coefficients` `(c₀, c₁, ..., cₙ)`, ordered from the lowest
to the highest degree, and its first derivative at `x`, both computed with the Horner
scheme. The tuple `coefficients` must have at least two elements.

# Returns

- `Number`: Polynomial value `c₀ + c₁ x + ... + cₙ xⁿ`.
- `Number`: Polynomial derivative `c₁ + 2 c₂ x + ... + n cₙ xⁿ⁻¹`.
"""
@inline function _evalpoly_with_derivative(
    x::Number, coefficients::NTuple{N, Number}
) where {N}
    p  = evalpoly(x, coefficients)
    ∂p = evalpoly(x, ntuple(i -> i * coefficients[i + 1], Val(N - 1)))
    return p, ∂p
end
