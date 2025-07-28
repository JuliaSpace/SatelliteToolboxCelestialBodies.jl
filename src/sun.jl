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

function sun_position_mod(jd_tdb::Number)
    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # Mean anomaly of the Sun [deg].
    Ms = 357.529_109_2 + 35_999.050_34t_tdb

    # Convert Ms to [rad] and limit to the interval [0,2π].
    Ms = mod2pi(deg2rad(Ms))

    # Compute auxiliary variables.
    sin_Ms,  cos_Ms  = sincos(Ms)

    # Let's compute `sin(2Ms)` and `cos(2Ms)` more efficiently.
    sin_2Ms = 2 * sin_Ms * cos_Ms
    cos_2Ms = cos_Ms * cos_Ms - sin_Ms * sin_Ms

    # Mean longitude of the Sun [deg].
    λ_m = 280.460 + 36_000.771t_tdb

    # Ecliptic latitude of the Sun [deg].
    λ_e = λ_m + 1.914_666_471sin_Ms + 0.019_994_643sin_2Ms

    # Obliquity of the ecliptic [deg].
    ϵ = 23.439_291 - 0.013_004_2t_tdb

    # Convert λ_e and ϵ to [rad] and limit to the interval [0,2π].
    λ_e = mod2pi(deg2rad(λ_e))
    ϵ   = mod2pi(deg2rad(ϵ))

    # Auxiliary variables.
    sin_ϵ   , cos_ϵ   = sincos(ϵ)
    sin_λ_e , cos_λ_e = sincos(λ_e)

    # Distance of the Sun from Earth [m].
    r = (1.000_140_612 - 0.016_708_617cos_Ms - 0.000_139_589cos_2Ms) * ASTRONOMICAL_UNIT

    # Compute the Sun vector represented in the Mean Equinox of Date (MOD).
    s_mod = SVector{3}(r * cos_λ_e, r * cos_ϵ * sin_λ_e, r * sin_ϵ * sin_λ_e)

    return s_mod
end

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

function sun_velocity_mod(jd_tdb::Number)
    # Convert centuries to seconds.
    cen2s = 1 / (36525.0 * 86400.0)

    # == Sun Position ======================================================================

    # TODO: Check if we can split this algorithm in a new function.
    #
    # This is the same algorithm as in `sun_position_mod`. However, we must have access to
    # some internal variables. Hence, we copied the entire function for now.

    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # Mean anomaly of the Sun [deg].
    Ms = 357.529_109_2 + 35_999.050_34t_tdb

    # Convert Ms to [rad] and limit to the interval [0,2π].
    Ms = mod2pi(deg2rad(Ms))

    # Compute auxiliary variables.
    sin_Ms,  cos_Ms  = sincos(Ms)

    # Let's compute `sin(2Ms)` and `cos(2Ms)` more efficiently.
    sin_2Ms = 2 * sin_Ms * cos_Ms
    cos_2Ms = cos_Ms * cos_Ms - sin_Ms * sin_Ms

    # Mean longitude of the Sun [deg].
    λ_m = 280.460 + 36_000.771t_tdb

    # Ecliptic latitude of the Sun [deg].
    λ_e = λ_m + 1.914_666_471sin_Ms + 0.019_994_643sin_2Ms

    # Obliquity of the ecliptic [deg].
    ϵ = 23.439_291 - 0.013_004_2t_tdb

    # Convert λ_e and ϵ to [rad] and limit to the interval [0,2π].
    λ_e = mod2pi(deg2rad(λ_e))
    ϵ   = mod2pi(deg2rad(ϵ))

    # Auxiliary variables.
    sin_ϵ   , cos_ϵ   = sincos(ϵ)
    sin_λ_e , cos_λ_e = sincos(λ_e)

    # Distance of the Sun from Earth [m].
    r = (1.000_140_612 - 0.016_708_617cos_Ms - 0.000_139_589cos_2Ms) * ASTRONOMICAL_UNIT

    # == Sun Velocity ======================================================================

    # Compute the required time derivatives [rad/s].
    ∂λ_m = deg2rad(36_000.771)   * cen2s
    ∂Ms  = deg2rad(35_999.05034) * cen2s
    ∂ϵ   = deg2rad(-0.0130042)  * cen2s
    ∂r   = (+0.016708617sin_Ms * ∂Ms + 2 * 0.000139589cos_2Ms * ∂Ms) * ASTRONOMICAL_UNIT
    ∂λ_e = ∂λ_m + deg2rad(1.914666471cos_Ms * ∂Ms + 2 * 0.019994643cos_2Ms * ∂Ms)

    # Compute the Sun velocity vector represented in the Mean Equinox of Date
    # (MOD).
    vsun_mod = SVector{3}(
        ∂r * cos_λ_e -                                    r * sin_λ_e * ∂λ_e,
        ∂r * cos_ϵ * sin_λ_e - r * sin_ϵ * sin_λ_e * ∂ϵ + r * cos_ϵ * cos_λ_e * ∂λ_e,
        ∂r * sin_ϵ * sin_λ_e + r * cos_ϵ * sin_λ_e * ∂ϵ + r * sin_ϵ * cos_λ_e * ∂λ_e
    )

    return vsun_mod
end
