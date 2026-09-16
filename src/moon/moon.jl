## Description #############################################################################
#
# Compute the Moon position and velocity.
#
## References ##############################################################################
#
# [1] Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
# [2] Meeus, J. (1998). Astronomical Algorithms. 2nd ed. Willmann-Bell, Inc, Richmond, VA.
#
############################################################################################

export moon_position_mod, moon_velocity_mod, moon_state_mod

############################################################################################
#                                      Moon Position                                       #
############################################################################################

"""
    moon_position_mod(jd_tdb::Number[, model]) -> SVector{3, T}
    moon_position_mod(date_tdb::Union{Date, DateTime}[, model]) -> SVector{3, Float64}

Compute the Moon position [m] represented in the IAU-76/FK5 MOD (mean-equator,
mean-equinox of date) at the Julian Day `jd_tdb` or at the date `date_tdb`, both in the
Barycentric Dynamical Time (TDB).

The `model` selects the algorithm and must be `Val(:Meeus)` **[2, ch. 47]** or
`Val(:Vallado)` **[1, p. 288]**. See [`moon_state_mod`](@ref) for the details on the
accuracy of each model and on the element type `T` of the result. The function throws an
error if `model` is not one of the supported models.

See also: [`moon_state_mod`](@ref), [`moon_velocity_mod`](@ref)

# Arguments

- `jd_tdb::Number`: Julian Day [TDB] at which the position must be computed.
- `date_tdb::Union{Date, DateTime}`: Date [TDB] at which the position must be computed.
- `model::Val`: Algorithm used to compute the Moon position.
    (**Default**: `Val(:Meeus)`)

# Returns

- `SVector{3, T}`: Moon position vector [m] represented in MOD.

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
- **[2]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.

# Extended help

## Throws

- `ArgumentError`: `model` is neither `Val(:Meeus)` nor `Val(:Vallado)`.
"""
function moon_position_mod(date_tdb::Union{Date, DateTime})
    return moon_position_mod(date_tdb, Val(:Meeus))
end

moon_position_mod(jd_tdb::Number) = moon_position_mod(jd_tdb, Val(:Meeus))

function moon_position_mod(date_tdb::Union{Date, DateTime}, model::Val)
    return moon_position_mod(datetime2julian(DateTime(date_tdb)), model)
end

# NOTE: Computing the velocity together with the position adds only a few operations per
# term of the series. Hence, we do not keep a separate position-only kernel.
moon_position_mod(jd_tdb::Number, model::Val) = moon_state_mod(jd_tdb, model)[1]

############################################################################################
#                                      Moon Velocity                                       #
############################################################################################

"""
    moon_velocity_mod(jd_tdb::Number[, model]) -> SVector{3, T}
    moon_velocity_mod(date_tdb::Union{Date, DateTime}[, model]) -> SVector{3, Float64}

Compute the Moon velocity [m/s] measured and represented in the IAU-76/FK5 MOD
(mean-equator, mean-equinox of date) at the Julian Day `jd_tdb` or at the date `date_tdb`,
both in the Barycentric Dynamical Time (TDB).

The `model` selects the algorithm and must be `Val(:Meeus)` **[2, ch. 47]** or
`Val(:Vallado)` **[1, p. 288]**. In both cases, the velocity is the analytical time
derivative of the Moon position. See [`moon_state_mod`](@ref) for the details on the
accuracy of each model and on the element type `T` of the result. The function throws an
error if `model` is not one of the supported models.

See also: [`moon_state_mod`](@ref), [`moon_position_mod`](@ref)

# Arguments

- `jd_tdb::Number`: Julian Day [TDB] at which the velocity must be computed.
- `date_tdb::Union{Date, DateTime}`: Date [TDB] at which the velocity must be computed.
- `model::Val`: Algorithm used to compute the Moon velocity.
    (**Default**: `Val(:Meeus)`)

# Returns

- `SVector{3, T}`: Moon velocity vector [m/s] represented in MOD.

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
- **[2]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.

# Extended help

## Throws

- `ArgumentError`: `model` is neither `Val(:Meeus)` nor `Val(:Vallado)`.
"""
function moon_velocity_mod(date_tdb::Union{Date, DateTime})
    return moon_velocity_mod(date_tdb, Val(:Meeus))
end

moon_velocity_mod(jd_tdb::Number) = moon_velocity_mod(jd_tdb, Val(:Meeus))

function moon_velocity_mod(date_tdb::Union{Date, DateTime}, model::Val)
    return moon_velocity_mod(datetime2julian(DateTime(date_tdb)), model)
end

moon_velocity_mod(jd_tdb::Number, model::Val) = moon_state_mod(jd_tdb, model)[2]

############################################################################################
#                                       Moon State                                        #
############################################################################################

"""
    moon_state_mod(jd_tdb::Number[, model]) -> SVector{3, T}, SVector{3, T}
    moon_state_mod(
        date_tdb::Union{Date, DateTime}[, model]
    ) -> SVector{3, Float64}, SVector{3, Float64}

Compute the Moon position [m] and velocity [m/s] represented in the IAU-76/FK5 MOD
(mean-equator, mean-equinox of date) at the Julian Day `jd_tdb` or at the date `date_tdb`,
both in the Barycentric Dynamical Time (TDB).

The `model` selects the algorithm and must be `Val(:Meeus)` or `Val(:Vallado)`.
`Val(:Meeus)` uses the algorithm in **[2, ch. 47]** that provides an accuracy of 10 [arcsec]
in the longitude and 4 [arcsec] in the latitude (the reference does not mention the
timespan). `Val(:Vallado)` uses the algorithm in **[1, p. 288]** that is about 3 times
faster than `Val(:Meeus)` but can lead to errors of 0.3 [°] in longitude and 0.2 [°] in
latitude. In both cases, the velocity is the analytical time derivative of the position.
Prefer this function over calling [`moon_position_mod`](@ref) and
[`moon_velocity_mod`](@ref) separately when both quantities are required, since the
computation is shared.

The function throws an error if `model` is not one of the supported models.

See also: [`moon_position_mod`](@ref), [`moon_velocity_mod`](@ref)

# Arguments

- `jd_tdb::Number`: Julian Day [TDB] at which the state must be computed.
- `date_tdb::Union{Date, DateTime}`: Date [TDB] at which the state must be computed.
- `model::Val`: Algorithm used to compute the Moon state.
    (**Default**: `Val(:Meeus)`)

# Returns

- `SVector{3, T}`: Moon position vector [m] represented in MOD.
- `SVector{3, T}`: Moon velocity vector [m/s] represented in MOD.

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
- **[2]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.

# Extended help

The element type `T` of the result is `promote_type(float(typeof(jd_tdb)), Float64)`.
Hence, `Float32` inputs yield `Float64` results because the precision of the series
requires it, whereas wider types, such as `BigFloat` or automatic differentiation numbers,
propagate to the output.

## Throws

- `ArgumentError`: `model` is neither `Val(:Meeus)` nor `Val(:Vallado)`.
"""
function moon_state_mod(date_tdb::Union{Date, DateTime})
    return moon_state_mod(date_tdb, Val(:Meeus))
end

moon_state_mod(jd_tdb::Number) = moon_state_mod(jd_tdb, Val(:Meeus))

function moon_state_mod(date_tdb::Union{Date, DateTime}, model::Val)
    return moon_state_mod(datetime2julian(DateTime(date_tdb)), model)
end

function moon_state_mod(::Number, ::Val{M}) where {M}
    return throw(
        ArgumentError(
            "The Moon model :$M is not supported. " *
            "The available models are :Meeus and :Vallado.",
        ),
    )
end

# NOTE: The formatter is disabled for the following functions to keep the hand-aligned
# numeric expressions, which improve the readability of the algorithms.
#! format: off
function moon_state_mod(jd_tdb::Number, ::Val{:Meeus})
    # Number of Julian centuries from the J2000 epoch [TDB].
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # == Fundamental Arguments =============================================================

    # The fundamental arguments are polynomials in `t_tdb` [deg], and their time derivatives
    # [deg/century] are obtained by differentiating the polynomials.

    # Moon's mean longitude referred to the mean equinox of date.
    L´_deg, ∂L´_deg = _evalpoly_with_derivative(
        t_tdb,
        (
            +218.316_447_7,
            +481_267.881_234_21,
            -0.001_578_6,
            +1 / 538_841,
            -1 / 65_194_000,
        )
    )

    # Mean elongation of the Moon.
    D_deg, ∂D_deg = _evalpoly_with_derivative(
        t_tdb,
        (
            +297.850_192_1,
            +445_267.111_403_4,
            -0.001_881_9,
            +1 / 545_868,
            -1 / 113_065_000,
        )
    )

    # Sun's mean anomaly.
    M_deg, ∂M_deg = _evalpoly_with_derivative(
        t_tdb,
        (
            +357.529_109_2,
            +35_999.050_290_9,
            -0.000_153_6,
            +1 / 24_490_000,
        )
    )

    # Moon's mean anomaly.
    M´_deg, ∂M´_deg = _evalpoly_with_derivative(
        t_tdb,
        (
            +134.963_396_4,
            +477_198.867_505_5,
            +0.008_741_4,
            +1 / 69_699,
            -1 / 14_712_000,
        )
    )

    # Moon's argument of latitude (mean distance of the Moon from its ascending node).
    F_deg, ∂F_deg = _evalpoly_with_derivative(
        t_tdb,
        (
            +93.272_095_0,
            +483_202.017_523_3,
            -0.003_653_9,
            -1 / 3_526_000,
            +1 / 863_310_000,
        )
    )

    # Obliquity of the ecliptic.
    ϵ_deg, ∂ϵ_deg = _evalpoly_with_derivative(
        t_tdb,
        (23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)
    )

    # Additional arguments required for the algorithm.
    A₁_deg, ∂A₁_deg = _evalpoly_with_derivative(t_tdb, (119.75, 131.849))
    A₂_deg, ∂A₂_deg = _evalpoly_with_derivative(t_tdb, ( 53.09, 479_264.290))
    A₃_deg, ∂A₃_deg = _evalpoly_with_derivative(t_tdb, (313.45, 481_266.484))

    # Convert the angles to [rad]. They are not wrapped to [0, 2π] because they are only
    # used as arguments of `sincos`, which performs its own range reduction.
    L´ = deg2rad(L´_deg)
    D  = deg2rad(D_deg)
    M  = deg2rad(M_deg)
    M´ = deg2rad(M´_deg)
    F  = deg2rad(F_deg)
    ϵ  = deg2rad(ϵ_deg)
    A₁ = deg2rad(A₁_deg)
    A₂ = deg2rad(A₂_deg)
    A₃ = deg2rad(A₃_deg)

    # Convert the angular rates to [rad/century].
    ∂L´ = deg2rad(∂L´_deg)
    ∂D  = deg2rad(∂D_deg)
    ∂M  = deg2rad(∂M_deg)
    ∂M´ = deg2rad(∂M´_deg)
    ∂F  = deg2rad(∂F_deg)
    ∂A₁ = deg2rad(∂A₁_deg)
    ∂A₂ = deg2rad(∂A₂_deg)
    ∂A₃ = deg2rad(∂A₃_deg)

    # Term used to correct the arguments that depend on the Sun's mean anomaly `M` due to
    # the decrease of the Earth's orbit eccentricity, and its time derivative [1/century].
    E, ∂E = _evalpoly_with_derivative(t_tdb, (1, -0.002_516, -0.000_007_4))

    # == Periodic Terms ====================================================================

    # Sines and cosines of the fundamental arguments. The arguments of the periodic terms
    # are integer combinations of them, whose sines and cosines are obtained with the angle
    # addition formulas instead of calling `sincos` for each term. Besides being faster,
    # this approach avoids rounding the combinations, which can reach thousands of radians
    # far from J2000.
    sc_D  = sincos(D)
    sc_M  = sincos(M)
    sc_M´ = sincos(M´)
    sc_F  = sincos(F)

    # Sum the periodic terms in the tables 47.A and 47.B [2] for the longitude, distance,
    # and latitude of the Moon, together with their time derivatives. The table 47.B has
    # only sine terms, so its cosine sums are discarded.
    Σl, Σr, ∂Σl, ∂Σr = _sum_periodic_terms(
        Val(:Table47A),
        sc_D, sc_M, sc_M´, sc_F, ∂D, ∂M, ∂M´, ∂F, E, ∂E
    )

    Σb, _, ∂Σb, _ = _sum_periodic_terms(
        Val(:Table47B),
        sc_D, sc_M, sc_M´, sc_F, ∂D, ∂M, ∂M´, ∂F, E, ∂E
    )

    # Apply the additive terms due to the action of Venus (`A₁`), Jupiter (`A₂`), and the
    # flattening of the Earth (`L´`) [2, p. 338], together with their time derivatives.
    sc_A₁ = sincos(A₁)
    sc_A₂ = sincos(A₂)
    sc_A₃ = sincos(A₃)
    sc_L´ = sincos(L´)

    sin_A₁,     cos_A₁     = sc_A₁
    sin_A₂,     cos_A₂     = sc_A₂
    sin_A₃,     cos_A₃     = sc_A₃
    sin_L´,     cos_L´     = sc_L´
    sin_L´_mF,  cos_L´_mF  = _sincos_sub(sc_L´, sc_F)
    sin_A₁_mF,  cos_A₁_mF  = _sincos_sub(sc_A₁, sc_F)
    sin_A₁_pF,  cos_A₁_pF  = _sincos_add(sc_A₁, sc_F)
    sin_L´_mM´, cos_L´_mM´ = _sincos_sub(sc_L´, sc_M´)
    sin_L´_pM´, cos_L´_pM´ = _sincos_add(sc_L´, sc_M´)

    Σl += 3958sin_A₁ + 1962sin_L´_mF + 318sin_A₂

    Σb += -2235sin_L´ +
        382sin_A₃ +
        175sin_A₁_mF +
        175sin_A₁_pF +
        127sin_L´_mM´ -
        115sin_L´_pM´

    ∂Σl += 3958cos_A₁ * ∂A₁ + 1962cos_L´_mF * (∂L´ - ∂F) + 318cos_A₂ * ∂A₂

    ∂Σb += -2235cos_L´ * ∂L´ +
        382cos_A₃ * ∂A₃ +
        175cos_A₁_mF * (∂A₁ - ∂F) +
        175cos_A₁_pF * (∂A₁ + ∂F) +
        127cos_L´_mM´ * (∂L´ - ∂M´) -
        115cos_L´_pM´ * (∂L´ + ∂M´)

    # == Moon Coordinates ==================================================================

    # Geocentric ecliptic longitude [rad] and latitude [rad] of the Moon referred to the
    # mean equinox of date, and its distance from the Earth's center [m]. The sums are
    # given in 10⁻⁶ deg and in m.
    λ = L´ + deg2rad(Σl / 1_000_000)
    β = deg2rad(Σb / 1_000_000)
    Δ = 385_000.56e3 + Σr

    # Time derivatives of the Moon coordinates [rad/s] and [m/s].
    ∂λ = (∂L´ + deg2rad(∂Σl / 1_000_000)) * _CENTURIES_PER_SECOND
    ∂β = deg2rad(∂Σb / 1_000_000) * _CENTURIES_PER_SECOND
    ∂Δ = ∂Σr * _CENTURIES_PER_SECOND
    ∂ϵ = deg2rad(∂ϵ_deg) * _CENTURIES_PER_SECOND

    # == Position and Velocity =============================================================

    # `λ` and `β` are the geocentric longitude and latitude of the Moon w.r.t. the mean
    # ecliptic and equinox of date. Hence, we must rotate the vector from the ecliptic to
    # the mean equator of date to obtain its representation in MOD.
    sin_λ, cos_λ = sincos(λ)
    sin_β, cos_β = sincos(β)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    # Moon position vector represented in MOD [m].
    r_moon_mod = SVector{3}(
        Δ * cos_β * cos_λ,
        Δ * (cos_ϵ * cos_β * sin_λ - sin_ϵ * sin_β),
        Δ * (sin_ϵ * cos_β * sin_λ + cos_ϵ * sin_β)
    )

    # Moon velocity vector represented in MOD [m/s], obtained by differentiating
    # `r_moon_mod`.
    v_moon_mod = SVector{3}(
        ∂Δ * cos_β * cos_λ -
            Δ * ∂β * sin_β * cos_λ -
            Δ * ∂λ * cos_β * sin_λ,
        ∂Δ * (cos_ϵ * cos_β * sin_λ - sin_ϵ * sin_β) -
            Δ * ∂ϵ * (sin_ϵ * cos_β * sin_λ + cos_ϵ * sin_β) -
            Δ * ∂β * (cos_ϵ * sin_β * sin_λ + sin_ϵ * cos_β) +
            Δ * ∂λ * cos_ϵ * cos_β * cos_λ,
        ∂Δ * (sin_ϵ * cos_β * sin_λ + cos_ϵ * sin_β) +
            Δ * ∂ϵ * (cos_ϵ * cos_β * sin_λ - sin_ϵ * sin_β) +
            Δ * ∂β * (cos_ϵ * cos_β - sin_ϵ * sin_β * sin_λ) +
            Δ * ∂λ * sin_ϵ * cos_β * cos_λ
    )

    return r_moon_mod, v_moon_mod
end

function moon_state_mod(jd_tdb::Number, ::Val{:Vallado})
    # Number of Julian centuries from the J2000 epoch [TDB].
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # == Fundamental Arguments =============================================================

    # Sines and cosines of the arguments of the periodic terms. The angles are not wrapped
    # to [0, 2π] because `sincos` performs its own range reduction.
    sin₁, cos₁ = sincos(deg2rad(134.9 + 477_198.85t_tdb))
    sin₂, cos₂ = sincos(deg2rad(259.2 - 413_335.38t_tdb))
    sin₃, cos₃ = sincos(deg2rad(235.7 + 890_534.23t_tdb))
    sin₄, cos₄ = sincos(deg2rad(269.9 + 954_397.70t_tdb))
    sin₅, cos₅ = sincos(deg2rad(357.5 +  35_999.05t_tdb))
    sin₆, cos₆ = sincos(deg2rad(186.6 + 966_404.05t_tdb))
    sin₇, cos₇ = sincos(deg2rad( 93.3 + 483_202.03t_tdb))
    sin₈, cos₈ = sincos(deg2rad(228.2 + 960_400.87t_tdb))
    sin₉, cos₉ = sincos(deg2rad(318.3 +   6_003.18t_tdb))
    sin₁₀, cos₁₀ = sincos(deg2rad(217.6 - 407_332.20t_tdb))

    # Ecliptic longitude of the Moon [deg].
    λₑ = 218.32 +
        481_267.8813t_tdb +
        6.29sin₁ -
        1.27sin₂ +
        0.66sin₃ +
        0.21sin₄ -
        0.19sin₅ -
        0.11sin₆

    # Ecliptic latitude of the Moon [deg].
    ϕₑ = 5.13sin₇ + 0.28sin₈ - 0.28sin₉ - 0.17sin₁₀

    # Horizontal parallax of the Moon [deg].
    P = 0.9508 + 0.0518cos₁ + 0.0095cos₂ + 0.0078cos₃ + 0.0028cos₄

    # Obliquity of the ecliptic [deg] and its time derivative [deg/century].
    ϵ, ∂ϵ = _evalpoly_with_derivative(
        t_tdb,
        (23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)
    )

    # Time derivative of the ecliptic longitude [deg/century].
    ∂λₑ = 481_267.8813 +
        6.29cos₁ * deg2rad( 477_198.85) -
        1.27cos₂ * deg2rad(-413_335.38) +
        0.66cos₃ * deg2rad( 890_534.23) +
        0.21cos₄ * deg2rad( 954_397.70) -
        0.19cos₅ * deg2rad(  35_999.05) -
        0.11cos₆ * deg2rad( 966_404.05)

    # Time derivative of the ecliptic latitude [deg/century].
    ∂ϕₑ = 5.13cos₇ * deg2rad( 483_202.03) +
        0.28cos₈ * deg2rad( 960_400.87) -
        0.28cos₉ * deg2rad(   6_003.18) -
        0.17cos₁₀ * deg2rad(-407_332.20)

    # Time derivative of the horizontal parallax [deg/century].
    ∂P = -0.0518sin₁ * deg2rad( 477_198.85) -
        0.0095sin₂ * deg2rad(-413_335.38) -
        0.0078sin₃ * deg2rad( 890_534.23) -
        0.0028sin₄ * deg2rad( 954_397.70)

    # Convert the angles to [rad] and the rates to [rad/s].
    λₑ = deg2rad(λₑ)
    ϕₑ = deg2rad(ϕₑ)
    P  = deg2rad(P)
    ϵ  = deg2rad(ϵ)

    ∂λₑ = deg2rad(∂λₑ) * _CENTURIES_PER_SECOND
    ∂ϕₑ = deg2rad(∂ϕₑ) * _CENTURIES_PER_SECOND
    ∂P  = deg2rad(∂P)  * _CENTURIES_PER_SECOND
    ∂ϵ  = deg2rad(∂ϵ)  * _CENTURIES_PER_SECOND

    # == Moon Coordinates ==================================================================

    sin_λ, cos_λ = sincos(λₑ)
    sin_ϕ, cos_ϕ = sincos(ϕₑ)
    sin_P, cos_P = sincos(P)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    # Distance from the Earth's center to the Moon [m], obtained from the horizontal
    # parallax, and its time derivative [m/s].
    r  = WGS84_ELLIPSOID.a / sin_P
    ∂r = -r * (cos_P / sin_P) * ∂P

    # == Position and Velocity =============================================================

    # `λₑ` and `ϕₑ` are the ecliptic longitude and latitude of the Moon w.r.t. the mean
    # equinox of date. Hence, we must rotate the vector from the ecliptic to the mean
    # equator of date to obtain its representation in MOD.

    # Moon position vector represented in MOD [m].
    r_moon_mod = SVector{3}(
        r * (cos_ϕ * cos_λ),
        r * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ),
        r * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ)
    )

    # Moon velocity vector represented in MOD [m/s], obtained by differentiating
    # `r_moon_mod`.
    v_moon_mod = SVector{3}(
        ∂r * cos_ϕ * cos_λ -
            r * ∂ϕₑ * sin_ϕ * cos_λ -
            r * ∂λₑ * cos_ϕ * sin_λ,
        ∂r * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ) -
            r * ∂ϵ * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ) -
            r * ∂ϕₑ * (cos_ϵ * sin_ϕ * sin_λ + sin_ϵ * cos_ϕ) +
            r * ∂λₑ * cos_ϵ * cos_ϕ * cos_λ,
        ∂r * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ) +
            r * ∂ϵ * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ) +
            r * ∂ϕₑ * (cos_ϵ * cos_ϕ - sin_ϵ * sin_ϕ * sin_λ) +
            r * ∂λₑ * sin_ϵ * cos_ϕ * cos_λ
    )

    return r_moon_mod, v_moon_mod
end
#! format: on

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _periodic_terms_table(::Val{:Table47A}) -> NTuple{60, SVector{6, Int32}}
    _periodic_terms_table(::Val{:Table47B}) -> NTuple{60, SVector{5, Int32}}

Return the table of periodic terms selected by the tag: [`_TAB_47A`](@ref) for
`Val(:Table47A)` and [`_TAB_47B`](@ref) for `Val(:Table47B)`.
"""
_periodic_terms_table(::Val{:Table47A}) = _TAB_47A
_periodic_terms_table(::Val{:Table47B}) = _TAB_47B

"""
    _sum_periodic_terms(
        ::Val{table},
        sc_D::Tuple{T, T},
        sc_M::Tuple{T, T},
        sc_M´::Tuple{T, T},
        sc_F::Tuple{T, T},
        ∂D::Number,
        ∂M::Number,
        ∂M´::Number,
        ∂F::Number,
        E::Number,
        ∂E::Number
    ) -> Number, Number, Number, Number

Sum the periodic terms of the `table`, which must be `:Table47A` **[2, pp. 339-340]** or
`:Table47B` **[2, p. 341]**, given the tuples `sc_D`, `sc_M`, `sc_M´`, and `sc_F` with the
sines and cosines of the fundamental arguments `D`, `M`, `M´`, and `F`, their time
derivatives `∂D`, `∂M`, `∂M´`, and `∂F` [rad/century], the eccentricity correction factor
`E` [-], and its time derivative `∂E` [1/century].

The sums are computed as:

```math
Σ_s = \\sum_k c_{s, k} E_k \\sin(θ_k), \\quad Σ_c = \\sum_k c_{c, k} E_k \\cos(θ_k),
```

where `θ_k` is the integer combination of the fundamental arguments of the term `k`,
`c_{s, k}` and `c_{c, k}` are its sine and cosine coefficients, and `E_k` is `E`, `E²`, or
1 depending on the multiplier of `M`. The tables without cosine coefficients yield zero
cosine sums.

The function is generated from the selected table so that the multipliers become literal
constants: the sines and cosines of the arguments are composed from those of the
fundamental arguments using the angle addition formulas, and the terms with null
multipliers are removed at compile time. Hence, the sines and cosines in `sc_D`, `sc_M`,
`sc_M´`, and `sc_F` must be accurate since every term is derived from them.

# Returns

- `Number`: Sum of the sine terms Σₛ, in the unit of the sine coefficients of the table.
- `Number`: Sum of the cosine terms Σ_c, in the unit of the cosine coefficients of the
    table.
- `Number`: Time derivative of Σₛ [1/century].
- `Number`: Time derivative of Σ_c [1/century].

# References

- **[2]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.
"""
@generated function _sum_periodic_terms(
    ::Val{table},
    sc_D::Tuple{T, T},
    sc_M::Tuple{T, T},
    sc_M´::Tuple{T, T},
    sc_F::Tuple{T, T},
    ∂D::Number,
    ∂M::Number,
    ∂M´::Number,
    ∂F::Number,
    E::Number,
    ∂E::Number,
) where {table, T <: Number}
    # NOTE: This function is generated because the table is only known at compile time
    # through its tag, and we want the multipliers to be literal constants in the
    # unrolled code, which lets the compiler remove the null terms and the runtime
    # selection of the eccentricity correction.
    rows = _periodic_terms_table(Val(table))

    # Maximum absolute multiplier of each fundamental argument, defining how many multiples
    # of its sine and cosine must be computed.
    n_D, n_M, n_M´, n_F = ntuple(i -> maximum(Int(abs(row[i])) for row in rows), 4)

    # Symbols of the local variables holding the multiples of each fundamental argument and
    # of the time derivatives of the arguments, in the same order of the table columns.
    multiples   = (:mult_D, :mult_M, :mult_M´, :mult_F)
    derivatives = (:∂D, :∂M, :∂M´, :∂F)

    terms = Expr[]

    for row in rows
        aD, aM, aM´, aF = Int(row[1]), Int(row[2]), Int(row[3]), Int(row[4])
        c_s = Int(row[5])
        c_c = length(row) == 6 ? Int(row[6]) : 0

        # Sine and cosine of the argument, composed only from the non-null multipliers.
        sc_factors = Expr[]
        ∂arg_terms = Expr[]

        for (a, mult, ∂) in zip((aD, aM, aM´, aF), multiples, derivatives)
            a == 0 && continue
            push!(sc_factors, :(_sincos_multiple($mult, $a)))
            push!(∂arg_terms, :($a * $∂))
        end

        sc_arg =
            isempty(sc_factors) ? :((zero(T), one(T))) :
            foldl((sc_a, sc_b) -> :(_sincos_add($sc_a, $sc_b)), sc_factors)

        ∂arg = isempty(∂arg_terms) ? :(zero(T)) : Expr(:call, :+, ∂arg_terms...)

        # Eccentricity correction factor of the term and its time derivative.
        E_k, ∂E_k = if abs(aM) == 1
            :E, :∂E
        elseif abs(aM) == 2
            :E², :∂E²
        else
            1, 0
        end

        push!(terms, quote
            sin_arg, cos_arg = $sc_arg
            ∂arg = $∂arg
            Σ_s += $c_s * $E_k * sin_arg
            ∂Σ_s += $c_s * ($∂E_k * sin_arg + $E_k * cos_arg * ∂arg)
        end)

        c_c == 0 && continue

        push!(terms, quote
            Σ_c  += $c_c * $E_k * cos_arg
            ∂Σ_c += $c_c * ($∂E_k * cos_arg - $E_k * sin_arg * ∂arg)
        end)
    end

    return quote
        R = promote_type(T, typeof(E))

        Σ_s  = zero(R)
        Σ_c  = zero(R)
        ∂Σ_s = zero(R)
        ∂Σ_c = zero(R)

        E²  = E * E
        ∂E² = 2E * ∂E

        # Sines and cosines of the multiples of the fundamental arguments.
        mult_D  = _sincos_multiples(sc_D, Val($n_D))
        mult_M  = _sincos_multiples(sc_M, Val($n_M))
        mult_M´ = _sincos_multiples(sc_M´, Val($n_M´))
        mult_F  = _sincos_multiples(sc_F, Val($n_F))

        $(terms...)

        return Σ_s, Σ_c, ∂Σ_s, ∂Σ_c
    end
end
