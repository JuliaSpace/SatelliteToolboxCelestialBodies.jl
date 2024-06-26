## Description #############################################################################
#
# Compute the Moon position.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
# [2] Meeus, J (1998). Astronomical algorithms. Willmann-Bell, Inc, Richmond, VA.
#
############################################################################################

export moon_position_mod

"""
    moon_position_mod(jd_tdb::Number[, model]) -> SVector{3, Float64}
    moon_position_mod(date_tdb::DateTime[, model]) -> SVector{3, Float64}

Compute the Moon position represented in the IAU-76/FK5 MOD (mean-equator, mean-equinox of
date) at the Julian Day `jd_tdb` or `date_tdb`. The input time must be represented in the
Barycentric Dynamical Time (TDB).

The `model` must be `Val(:Meeus)` or `Val(:Vallado)`. `Val(:Meeus)` uses the algorithm in
**[2, p. 337]** that provides an accuracy of 10" in the longitude and 4" in the latitude
(the reference does not mention the timespan). `Val(:Vallado)` uses the algorithm in
**[1, p. 288]** that is 10x faster than `Val(:Meeus)` but can lead to errors of 0.3° in
longitude and 0.2° in latitude.

!!! note

    This function performs all the computations using `Float64` due to the necessary
    precision.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
- **[2]** Meeus, J (1998). Astronomical algorithms. Willmann-Bell, Inc, Richmond, VA.
"""
moon_position_mod(date_tdb::DateTime) = moon_position_mod(date_tdb, Val(:Meeus))
moon_position_mod(jd_tdb::Number) = moon_position_mod(jd_tdb, Val(:Meeus))

function moon_position_mod(date_tdb::DateTime, ::Val{:Meeus})
    return moon_position_mod(datetime2julian(date_tdb), Val(:Meeus))
end

function moon_position_mod(jd_tdb::Number, ::Val{:Meeus})
    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # Moon's mean latitude referred to the mean equinox of data [deg].
    L´ = @evalpoly(
        t_tdb,
        +218.316_447_7,
        +481_267.881_234_21,
        -0.001_578_6,
        +1 / 538_841,
        -1 / 65_194_000
    )

    # Mean elongation of the Moon [deg].
    D = @evalpoly(
        t_tdb,
        +297.850_192_1,
        +445_267.111_403_4,
        -0.001_881_9,
        +1 / 545_868,
        -1 / 113_065_000
    )

    # Sun's mean anomaly [deg].
    M = @evalpoly(
        t_tdb,
        +357.529_109_2,
        +35_999.050_290_9,
        -0.000_153_6,
        +1 / 24_490_000
    )

    # Moon's mean anomaly [deg].
    M´ = @evalpoly(
        t_tdb,
        +134.963_396_4,
        +477_198.867_505_5,
        +0.008_741_4,
        +1 / 69_699,
        -1 / 14_712_000
    )

    # Moon's argument of latitude (mean distance of the Moon from its ascending node) [deg].
    F = @evalpoly(
        t_tdb,
        +93.272_095_0,
        +483_202.017_523_3,
        -0.003_653_9,
        -1 / 3_526_000,
        +1 / 863_310_000
    )

    # Obliquity of the ecliptic [deg].
    ϵ = @evalpoly(t_tdb, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)

    # Additional arguments required for the algorithm [deg].
    A₁ = @evalpoly(t_tdb, 119.75, 131.849)
    A₂ = @evalpoly(t_tdb,  53.09, 479_264.290)
    A₃ = @evalpoly(t_tdb, 313.45, 481_266.484)

    # Convert everything to [rad] and limit the angles between [0, 2π].
    L´ = mod2pi(deg2rad(L´))
    D  = mod2pi(deg2rad(D))
    M  = mod2pi(deg2rad(M))
    M´ = mod2pi(deg2rad(M´))
    F  = mod2pi(deg2rad(F))
    ϵ  = mod2pi(deg2rad(ϵ))
    A₁ = mod2pi(deg2rad(A₁))
    A₂ = mod2pi(deg2rad(A₂))
    A₃ = mod2pi(deg2rad(A₃))

    # Term used to correct the arguments of the angle M that depends on the Earth's orbit
    # eccentricity around the Sun.
    E  = @evalpoly(t_tdb, 1, -0.002516, 0.0000074)
    E² = E * E

    # Compute the sum of the terms in the tables 47.A and 47.B [2].

    # Auxiliary variables to simplify the expressions.
    tab = _TAB_47A

    # Longitude and distance of the Moon.
    num_terms = size(tab)[2]

    Σl = 0.0
    Σr = 0.0

    @inbounds for k in 1:num_terms
        aD  = tab[1, k]
        aM  = tab[2, k]
        aM´ = tab[3, k]
        aF  = tab[4, k]

        arg = mod2pi(aD * D + aM * M + aM´ * M´ + aF * F)

        # Check if we need to apply the correction `E`.
        E_corr = if ((aM == 1) || (aM == -1))
            E
        elseif ((aM == 2) || (aM == -2))
            E²
        else
            one(E)
        end

        # Compute the quantities.
        sin_arg, cos_arg = sincos(arg)

        Σl += tab[5, k] * E_corr * sin_arg
        Σr += tab[6, k] * E_corr * cos_arg
    end

    # Auxiliary variables to simplify the expressions.
    tab = _TAB_47B

    # Latitude of the Moon.
    num_terms = size(tab)[2]

    Σb = 0.0

    @inbounds for k in 1:num_terms
        aD  = tab[1, k]
        aM  = tab[2, k]
        aM´ = tab[3, k]
        aF  = tab[4, k]

        arg = mod2pi(aD * D + aM * M + aM´ * M´ + aF * F)

        # Check if we need to apply the correction `E`.
        E_corr = if ((aM == 1) || (aM == -1))
            E
        elseif ((aM == 2) || (aM == -2))
            E²
        else
            one(E)
        end

        Σb += tab[5, k] * E_corr * sin(arg)
    end

    # Apply the corrections to the terms.
    Σl +=  3958sin(A₁) + 1962sin(L´ - F) + 318sin(A₂)
    Σb += -2235sin(L´) + 382sin(A₃) + 175sin(A₁ - F) + 175sin(A₁ + F) + 127sin(L´ - M´) - 115sin(L´ + M´)

    # Convert to [rad].
    Σl = mod2pi(deg2rad(Σl / 1_000_000))
    Σb = mod2pi(deg2rad(Σb / 1_000_000))

    # Compute the Moon coordinates [rad] and [m].
    λ = mod2pi(L´ + Σl)
    β = Σb
    Δ = 385_000.56e3 + Σr

    # Compute the Moon vector in MOD [m].
    #
    # Notice that λ and β provide us the geocentric latitude and longitude of the Moon
    # w.r.t. the mean equator of date in the ecliptic plane. Hence, we need also to rotate
    # the mean ecliptic to obtain the vector in the MOD.
    sin_λ, cos_λ = sincos(λ)
    sin_β, cos_β = sincos(β)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    r_moon_mod = @SVector [
        Δ * cos_β * cos_λ,
        Δ * (cos_ϵ * cos_β * sin_λ - sin_ϵ * sin_β),
        Δ * (sin_ϵ * cos_β * sin_λ + cos_ϵ * sin_β)
    ]

    return r_moon_mod
end

function moon_position_mod(date_tdb::DateTime, ::Val{:Vallado})
    return moon_position_mod(datetime2julian(date_tdb), Val(:Vallado))
end

function moon_position_mod(jd_tdb::Number, ::Val{:Vallado})
    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000)/36525.0

    # Auxiliary computation to improve performance.
    sin1, cos1 = sincos(deg2rad(134.9 + 477_198.85t_tdb))
    sin2, cos2 = sincos(deg2rad(259.2 - 413_335.38t_tdb))
    sin3, cos3 = sincos(deg2rad(235.7 + 890_534.23t_tdb))
    sin4, cos4 = sincos(deg2rad(269.9 + 954_397.70t_tdb))
    sin5       = sin(deg2rad(357.5 +  35_999.05t_tdb))
    sin6       = sin(deg2rad(186.6 + 966_404.05t_tdb))

    # Ecliptic latitude of the Moon [deg].
    λₑ = 218.32 + 481_267.8813t_tdb + 6.29sin1 - 1.27sin2 + 0.66sin3 + 0.21sin4 - 0.19sin5 - 0.11sin6

    # Ecliptic longitude of the Moon [deg].
    ϕₑ = 5.13sin(deg2rad( 93.3 + 483_202.03t_tdb)) +
         0.28sin(deg2rad(228.2 + 960_400.87t_tdb)) -
         0.28sin(deg2rad(318.3 +   6_003.18t_tdb)) -
         0.17sin(deg2rad(217.6 - 407_332.20t_tdb))

    # Parallax [deg].
    P = 0.9508 + 0.0518cos1 + 0.0095cos2 + 0.0078cos3 + 0.0028cos4

    # Obliquity of the ecliptic [deg].
    ϵ = @evalpoly(t_tdb, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)

    # Convert to radians and limit to the interval [0,2π].
    λₑ = mod2pi(deg2rad(λₑ))
    ϕₑ = mod2pi(deg2rad(ϕₑ))
    P  = mod2pi(deg2rad(P))
    ϵ  = mod2pi(deg2rad(ϵ))

    # Compute the distance from Earth to the Moon [m].
    r = WGS84_ELLIPSOID.a / sin(P)

    # Auxiliary variables.
    sin_λ, cos_λ = sincos(λₑ)
    sin_ϕ, cos_ϕ = sincos(ϕₑ)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    # Compute the Moon vector represented in MOD (IAU-76/KF5 mean-equator,
    # mean-equinox of date).
    r_moon_mod = @SVector [
        r * (cos_ϕ * cos_λ),
        r * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ),
        r * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ)
    ]

    return r_moon_mod
end
