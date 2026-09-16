## Description #############################################################################
#
# Compute the Sun position and velocity.
#
## References ##############################################################################
#
# [1] Vallado, D. A. (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
############################################################################################

export sun_position_mod, sun_velocity_mod, sun_state_mod

"""
    sun_position_mod(jd_tdb::Number) -> SVector{3, T}
    sun_position_mod(date_tdb::Union{Date, DateTime}) -> SVector{3, Float64}

Compute the Sun position [m] represented in the IAU-76/FK5 MOD (mean-equator, mean-equinox
of date) at the Julian Day `jd_tdb` or at the date `date_tdb`, both in the Barycentric
Dynamical Time (TDB).

The algorithm is described in **[1, pp. 277-279]**. See [`sun_state_mod`](@ref) for the
details on the accuracy and on the element type `T` of the result.

See also: [`sun_state_mod`](@ref), [`sun_velocity_mod`](@ref)

# Returns

- `SVector{3, T}`: Sun position vector [m] represented in MOD.

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function sun_position_mod(date_tdb::Union{Date, DateTime})
    return sun_position_mod(datetime2julian(DateTime(date_tdb)))
end

# NOTE: Computing the velocity together with the position costs only a handful of
# floating-point operations. Hence, we do not keep a separate position-only kernel.
sun_position_mod(jd_tdb::Number) = sun_state_mod(jd_tdb)[1]

"""
    sun_velocity_mod(jd_tdb::Number) -> SVector{3, T}
    sun_velocity_mod(date_tdb::Union{Date, DateTime}) -> SVector{3, Float64}

Compute the Sun velocity [m/s] measured and represented in the IAU-76/FK5 MOD
(mean-equator, mean-equinox of date) at the Julian Day `jd_tdb` or at the date `date_tdb`,
both in the Barycentric Dynamical Time (TDB).

The velocity is the analytical time derivative of the Sun position described in
**[1, pp. 277-279]**. See [`sun_state_mod`](@ref) for the details on the accuracy and on
the element type `T` of the result.

See also: [`sun_state_mod`](@ref), [`sun_position_mod`](@ref)

# Returns

- `SVector{3, T}`: Sun velocity vector [m/s] represented in MOD.

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function sun_velocity_mod(date_tdb::Union{Date, DateTime})
    return sun_velocity_mod(datetime2julian(DateTime(date_tdb)))
end

sun_velocity_mod(jd_tdb::Number) = sun_state_mod(jd_tdb)[2]

"""
    sun_state_mod(jd_tdb::Number) -> SVector{3, T}, SVector{3, T}
    sun_state_mod(date_tdb::Union{Date, DateTime}) -> SVector{3, Float64}, SVector{3, Float64}

Compute the Sun position [m] and velocity [m/s] represented in the IAU-76/FK5 MOD
(mean-equator, mean-equinox of date) at the Julian Day `jd_tdb` or at the date `date_tdb`,
both in the Barycentric Dynamical Time (TDB).

The position follows the algorithm in **[1, pp. 277-279]** and the velocity is its
analytical time derivative. Prefer this function over calling [`sun_position_mod`](@ref)
and [`sun_velocity_mod`](@ref) separately when both quantities are required, since the
computation is shared.

See also: [`sun_position_mod`](@ref), [`sun_velocity_mod`](@ref)

# Returns

- `SVector{3, T}`: Sun position vector [m] represented in MOD.
- `SVector{3, T}`: Sun velocity vector [m/s] represented in MOD.

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.

# Extended help

The element type `T` of the result is `promote_type(float(typeof(jd_tdb)), Float64)`.
Hence, `Float32` inputs yield `Float64` results because the precision of the series
requires it, whereas wider types, such as `BigFloat` or automatic differentiation numbers,
propagate to the output.

The algorithm uses the number of Julian centuries in TDB for all fundamental arguments,
including the mean longitude of the Sun, which is defined in UT1 in **[1]**. The resulting
error is below 1'' and negligible for the accuracy of this model.
"""
function sun_state_mod(date_tdb::Union{Date, DateTime})
    return sun_state_mod(datetime2julian(DateTime(date_tdb)))
end

# NOTE: The formatter is disabled for the following function to keep the hand-aligned
# numeric expressions, which improve the readability of the algorithm.
#! format: off
function sun_state_mod(jd_tdb::Number)
    # Number of Julian centuries from the J2000 epoch [TDB].
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # == Sun Position ======================================================================

    # Mean anomaly of the Sun [deg].
    Ms = 357.529_109_2 + 35_999.050_34t_tdb

    # Convert `Ms` to [rad]. The angles are not wrapped to [0, 2π] because they are only
    # used as arguments of `sincos`, which performs its own range reduction.
    Ms = deg2rad(Ms)

    sin_Ms, cos_Ms = sincos(Ms)

    # Compute `sin(2Ms)` and `cos(2Ms)` using the double-angle formulas.
    sin_2Ms = 2 * sin_Ms * cos_Ms
    cos_2Ms = cos_Ms * cos_Ms - sin_Ms * sin_Ms

    # Mean longitude of the Sun [deg].
    λ_m = 280.460 + 36_000.771t_tdb

    # Ecliptic longitude of the Sun [deg].
    λ_e = λ_m + 1.914_666_471sin_Ms + 0.019_994_643sin_2Ms

    # Obliquity of the ecliptic [deg].
    ϵ = 23.439_291 - 0.013_004_2t_tdb

    # Convert `λ_e` and `ϵ` to [rad].
    λ_e = deg2rad(λ_e)
    ϵ   = deg2rad(ϵ)

    sin_ϵ,   cos_ϵ   = sincos(ϵ)
    sin_λ_e, cos_λ_e = sincos(λ_e)

    # Distance of the Sun from Earth [m].
    r = (1.000_140_612 - 0.016_708_617cos_Ms - 0.000_139_589cos_2Ms) * ASTRONOMICAL_UNIT

    # Sun position vector represented in MOD [m].
    s_mod = SVector{3}(r * cos_λ_e, r * cos_ϵ * sin_λ_e, r * sin_ϵ * sin_λ_e)

    # == Sun Velocity ======================================================================

    # Time derivatives of the fundamental arguments [rad/s].
    ∂λ_m = deg2rad(36_000.771)    * _CENTURIES_PER_SECOND
    ∂Ms  = deg2rad(35_999.050_34) * _CENTURIES_PER_SECOND
    ∂ϵ   = deg2rad(-0.013_004_2)  * _CENTURIES_PER_SECOND

    # Time derivatives of the Sun distance [m/s] and of the ecliptic longitude [rad/s].
    ∂r   = (0.016_708_617sin_Ms + 2 * 0.000_139_589sin_2Ms) * ∂Ms * ASTRONOMICAL_UNIT
    ∂λ_e = ∂λ_m + deg2rad(1.914_666_471cos_Ms + 2 * 0.019_994_643cos_2Ms) * ∂Ms

    # Sun velocity vector represented in MOD [m/s], obtained by differentiating `s_mod`.
    ṡ_mod = SVector{3}(
        ∂r * cos_λ_e - r * sin_λ_e * ∂λ_e,
        ∂r * cos_ϵ * sin_λ_e - r * sin_ϵ * sin_λ_e * ∂ϵ +
            r * cos_ϵ * cos_λ_e * ∂λ_e,
        ∂r * sin_ϵ * sin_λ_e + r * cos_ϵ * sin_λ_e * ∂ϵ +
            r * sin_ϵ * cos_λ_e * ∂λ_e
    )

    return s_mod, ṡ_mod
end
#! format: on
