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

"""
    _sincos_add(sc_a::Tuple{Number, Number}, sc_b::Tuple{Number, Number}) -> Number, Number

Compute the sine and cosine of the sum of two angles given the tuples `sc_a` and `sc_b`
containing the sine and cosine of each angle, using the angle addition formulas.

# Returns

- `Number`: Sine of the sum of the angles.
- `Number`: Cosine of the sum of the angles.
"""
@inline function _sincos_add(sc_a::Tuple{Number, Number}, sc_b::Tuple{Number, Number})
    sin_a, cos_a = sc_a
    sin_b, cos_b = sc_b
    return sin_a * cos_b + cos_a * sin_b, cos_a * cos_b - sin_a * sin_b
end

"""
    _sincos_sub(sc_a::Tuple{Number, Number}, sc_b::Tuple{Number, Number}) -> Number, Number

Compute the sine and cosine of the difference of two angles given the tuples `sc_a` and
`sc_b` containing the sine and cosine of the minuend and subtrahend angles, using the angle
subtraction formulas.

# Returns

- `Number`: Sine of the difference of the angles.
- `Number`: Cosine of the difference of the angles.
"""
@inline function _sincos_sub(sc_a::Tuple{Number, Number}, sc_b::Tuple{Number, Number})
    sin_b, cos_b = sc_b
    return _sincos_add(sc_a, (-sin_b, cos_b))
end

"""
    _sincos_multiples(sc_θ::Tuple{T, T}, ::Val{N}) -> NTuple{N + 1, Tuple{T, T}}

Compute the sine and cosine of the integer multiples `0θ`, `1θ`, ..., `Nθ` of an angle
given the tuple `sc_θ` with its sine and cosine, using the angle addition formulas.

The result is a tuple in which the element `k + 1` contains the sine and cosine of `kθ`.

# Returns

- `NTuple{N + 1, Tuple{T, T}}`: Sines and cosines of the multiples of the angle.
"""
@inline function _sincos_multiples(sc_θ::Tuple{T, T}, ::Val{N}) where {T <: Number, N}
    sc_0 = (zero(T), one(T))

    # NOTE: The loops are unrolled at compile time since `N` is small. The recurrence
    # accumulates one rounding error per step, which is negligible for the small
    # multipliers used by the algorithms.
    return ntuple(Val(N + 1)) do i
        # The element `i` holds the multiple `k = i - 1`, obtained by adding `θ` `k` times.
        sc_kθ = sc_0

        for _ in 1:(i - 1)
            sc_kθ = _sincos_add(sc_kθ, sc_θ)
        end

        sc_kθ
    end
end

"""
    _sincos_multiple(sc_multiples::Tuple{Tuple{T, T}, Vararg{Tuple{T, T}}}, k::Integer) -> T, T

Return the sine and cosine of `kθ` given the tuple `sc_multiples` with the sines and
cosines of the non-negative multiples of `θ`, as computed by [`_sincos_multiples`](@ref).
The multiplier `k` can be negative, in which case the parity of the sine and cosine
functions is used. `abs(k)` must be lower than the length of `sc_multiples`.

# Returns

- `T`: Sine of `kθ`.
- `T`: Cosine of `kθ`.
"""
@inline function _sincos_multiple(
    sc_multiples::Tuple{Tuple{T, T}, Vararg{Tuple{T, T}}}, k::Integer
) where {T <: Number}
    sin_kθ, cos_kθ = sc_multiples[abs(k) + 1]
    return (k < 0 ? -sin_kθ : sin_kθ), cos_kθ
end
