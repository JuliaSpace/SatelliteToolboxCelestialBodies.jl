# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Compute the Sun position and velocity.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] The Astronomical Almanac for the year 2000 (p. C24).
#
#   [2] http://aa.usno.navy.mil/faq/docs/SunApprox.php
#
#   [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#       Microcosm Press, Hawthorne, CA.
#
#   [4] The Astronomical Almanac for the year 2006.
#
#   [5] Blanco, Manuel Jesus, Milidonis, Kypros, Bonanos, Aristides. Updating the PSA sun
#       position algorithm. Solar Energy, vol.212, Elsevier BV,2020-12.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export sun_position_mod, sun_velocity_mod, sun_position_el

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

    # Mean long of the Sun [deg].
    λ_m = 280.460 + 36_000.771t_tdb

    # Ecliptic lat of the Sun [deg].
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
    s_mod = @SVector [r * cos_λ_e, r * cos_ϵ * sin_λ_e, r * sin_ϵ * sin_λ_e]

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

    # Sun position
    # ======================================================================================

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

    # Mean long of the Sun [deg].
    λ_m = 280.460 + 36_000.771t_tdb

    # Ecliptic lat of the Sun [deg].
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

    # Sun velocity
    # ==========================================================================

    # Compute the required time derivatives [rad/s].
    ∂λ_m = deg2rad(36_000.771)   * cen2s
    ∂Ms  = deg2rad(35_999.05034) * cen2s
    ∂ϵ   = deg2rad(-0.0130042)  * cen2s
    ∂r   = (+0.016708617sin_Ms * ∂Ms + 2 * 0.000139589cos_2Ms * ∂Ms) * ASTRONOMICAL_UNIT
    ∂λ_e = ∂λ_m + deg2rad(1.914666471cos_Ms * ∂Ms + 2 * 0.019994643cos_2Ms * ∂Ms)

    # Compute the Sun velocity vector represented in the Mean Equinox of Date
    # (MOD).
    vsun_mod = @SVector [
        ∂r * cos_λ_e -                                    r * sin_λ_e * ∂λ_e,
        ∂r * cos_ϵ * sin_λ_e - r * sin_ϵ * sin_λ_e * ∂ϵ + r * cos_ϵ * cos_λ_e * ∂λ_e,
        ∂r * sin_ϵ * sin_λ_e + r * cos_ϵ * sin_λ_e * ∂ϵ + r * sin_ϵ * cos_λ_e * ∂λ_e
    ]

    return vsun_mod
end

"""
    sun_position_el(
        jd::Number,
        lat::Real=0.0,
        long::Real=0.0,
        flag::Char='l',
    )

Compute the Sun position represented in the Local Horizon Reference System at the
Julian Day `jd`, lat `lat`, long `long`, Atmospheric pressure `Pressure`,
Ambient Temperature `Temperature`, Output Flag `flag` and Algorithm `algorithm`. It utilises
the PSA+ algorithm, (\\[5] (accessed on 2023-07-09)), to compute the topocentric sun position.

# Inputs

    - jd: Julian Day;
    - lat: Latitude of the observer, in degrees, WSG84;
    - long: Longitude of the observer, in degrees, WSG84;

# Outputs

## Equatorial system:
    - α, Right Ascension in degrees;
    - δ, Declination in degrees;

## Local Coordinates:
    - ω, Hour angle in degrees;
    - θ, Zenith in degrees;
    - γ, Azimuth in degrees;

## Sun vector in (East, North, Zenith):
    - SunVec, [Nx3] sun vector in (east, north, zenith);

# References

- **[5]**: Blanco, Manuel Jesus, Milidonis, Kypros, Bonanos, Aristides. Updating the PSA sun
            position algorithm. Solar Energy, vol.212, Elsevier BV,2020-12.

"""
function sun_position_el(
    jd::Real,
    lat::Real=0.0,
    long::Real=0.0,
)
    # Get time data from Julian Date `jd`
    elapsedJD = jd - JD_J2000
    _, _, _, Hour, Minute, Second = jd_to_date(jd)
    DecimalHours = Hour + Minute/60.0 + Second/3600.0

    # PSA+ Algorithm
    ## Ecliptic Coordinates
    Ω =  2.267127827e+00 - 9.300339267e-04*elapsedJD #
    ML = 4.895036035e+00 + 1.720279602e-02*elapsedJD # Mean long
    MA = 6.239468336e+00 + 1.720200135e-02*elapsedJD # Mean Anomaly

    # Ecliptic long
    λ₀ = (
        ML + 3.338320972e-02*sin( MA )
        + 3.497596876e-04 * sin( 2*MA ) - 1.544353226e-04
        - 8.689729360e-06*sin( Ω )
    )
    # Ecliptic Obliquity
    ϵ₀ = 4.090904909e-01 - 6.213605399e-09*elapsedJD + 4.418094944e-05*cos(Ω)

    ## Celestial coordinates
    # Right ascension & declination
    dY1 = cos( ϵ₀ ) .* sin( λ₀ )
    dX1 = cos( λ₀ )
    α = mod2pi(atan( dY1, dX1))
    δ = asin( sin( ϵ₀ ) .* sin( λ₀ ) )

    ## Topocentric coordinates
    # Greenwich & Local sidereal time
    GMST = 6.697096103e+00 + 6.570984737e-02*elapsedJD + DecimalHours
    LMST = deg2rad( GMST*15 + long )
    # Hour angle
    ω = mod2pi(LMST - α)

    # Local coordinates
    # Zenith
    θ = acos(cosd(lat)*cos( ω ).*cos( δ ) + sin( δ )*sind(lat))
    dY2 = -sin(ω)
    dX2 = tan( δ) * cosd( lat ) - sind(lat)*cos(ω)
    # Azimuth
    γ = mod2pi(atan(dY2, dX2))

    # Parallax correction
    # ToDo: Add radius of earth in Base constants? Then import and replace it below?
    θ += 6371.01e+03/ASTRONOMICAL_UNIT * sin(θ)

    # East North Zenith Frame
    SunVec = [
        sin(γ).*sin(θ),
        cos(γ).*sin(θ),
        cos(θ)
    ]

    return (rad2deg(α), rad2deg(δ), rad2deg(ω), rad2deg(θ), rad2deg(γ), SunVec)
end
