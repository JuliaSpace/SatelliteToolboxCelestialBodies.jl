## Description #############################################################################
#
# Compute the Sun position and velocity.
#
## References ##############################################################################
#
# [1] The Astronomical Almanac for the year 2000 (p. C24).
#
# [2] http://aa.usno.navy.mil/faq/docs/SunApprox.php
#
# [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
# [4] The Astronomical Almanac for the year 2006.
#
############################################################################################

export sun_position_mod, sun_velocity_mod

"""
    sun_position_mod(jd_tdb::Number) -> SVector{3, Float64}
    sun_position_mod(date_tdb::DateTime) -> SVector{3, Float64}

Compute the Sun position represented in the IAU-76/FK5 MOD (mean-equator, mean-equinox of
date) at the Julian Day `jd_tdb` or `date_tdb`. The input time must be represented in the
Barycentric Dynamical Time (TDB). The algorithm was adapted from [1, pp. 277-279].

!!! Note
    This function performs all the computations using `Float64` due to the necessary
    precision.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
sun_position_mod(date_tdb::DateTime) = sun_position_mod(datetime2julian(date_tdb))

# NOTE: Computing the velocity together with the position costs only a handful of
# floating-point operations. Hence, we do not keep a separate position-only kernel.
sun_position_mod(jd_tdb::Number) = _sun_state_mod(jd_tdb)[1]

"""
    sun_velocity_mod(jd_tdb::Number) -> SVector{3, Float64}
    sun_velocity_mod(date_tdb::DateTime) -> SVector{3, Float64}

Compute the Sun velocity measured and represented in the IAU-76/FK5 MOD (mean-equator,
mean-equinox of date) at the Julian Day `jd_tdb` or `date_tdb`. The input time must be
represented in the Barycentric Dynamical Time (TDB). The algorithm was obtained by computing
the time derivative of the Sun position in [1, p. 277-279].

!!! note
    This function performs all the computations using `Float64` due to the necessary
    precision.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
sun_velocity_mod(date_tdb::DateTime) = sun_velocity_mod(datetime2julian(date_tdb))
sun_velocity_mod(jd_tdb::Number) = _sun_state_mod(jd_tdb)[2]

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _sun_state_mod(jd_tdb::Number) -> SVector{3, Float64}, SVector{3, Float64}

Compute the Sun position [m] and velocity [m/s] represented in the IAU-76/FK5 MOD
(mean-equator, mean-equinox of date) at the Julian Day `jd_tdb` [TDB].

The position follows the algorithm in **[1, pp. 277-279]**. The velocity is its analytical
time derivative.

# Returns

- `SVector{3, Float64}`: Sun position vector [m] represented in MOD.
- `SVector{3, Float64}`: Sun velocity vector [m/s] represented in MOD.

# References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
function _sun_state_mod(jd_tdb::Number)
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
    ∂λ_m = deg2rad(36_000.771)    * _CENTURY_TO_SECONDS
    ∂Ms  = deg2rad(35_999.050_34) * _CENTURY_TO_SECONDS
    ∂ϵ   = deg2rad(-0.013_004_2)  * _CENTURY_TO_SECONDS

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
