## Description #############################################################################
#
# Compute the Moon position and velocity.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
# [2] Meeus, J (1998). Astronomical algorithms. Willmann-Bell, Inc, Richmond, VA.
#
############################################################################################

export moon_position_mod, moon_velocity_mod

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

    r_moon_mod = SVector{3}(
        Δ * cos_β * cos_λ,
        Δ * (cos_ϵ * cos_β * sin_λ - sin_ϵ * sin_β),
        Δ * (sin_ϵ * cos_β * sin_λ + cos_ϵ * sin_β)
    )

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
    r_moon_mod = SVector{3}(
        r * (cos_ϕ * cos_λ),
        r * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ),
        r * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ)
    )

    return r_moon_mod
end

############################################################################################
#                                      Moon Velocity                                       #
############################################################################################

"""
    moon_velocity_mod(jd_tdb::Number[, model]) -> SVector{3, Float64}
    moon_velocity_mod(date_tdb::DateTime[, model]) -> SVector{3, Float64}

Compute the Moon velocity measured and represented in the IAU-76/FK5 MOD (mean-equator,
mean-equinox of date) at the Julian Day `jd_tdb` or `date_tdb`. The input time must be
represented in the Barycentric Dynamical Time (TDB). The algorithm was obtained by computing
the time derivative of the Moon position.

The `model` must be `Val(:Meeus)` or `Val(:Vallado)`. See [`moon_position_mod`](@ref) for
details on the accuracy of each model.

!!! note

    This function performs all the computations using `Float64` due to the necessary
    precision.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
- **[2]** Meeus, J (1998). Astronomical algorithms. Willmann-Bell, Inc, Richmond, VA.
"""
moon_velocity_mod(date_tdb::DateTime) = moon_velocity_mod(date_tdb, Val(:Meeus))
moon_velocity_mod(jd_tdb::Number) = moon_velocity_mod(jd_tdb, Val(:Meeus))

function moon_velocity_mod(date_tdb::DateTime, ::Val{:Meeus})
    return moon_velocity_mod(datetime2julian(date_tdb), Val(:Meeus))
end

function moon_velocity_mod(jd_tdb::Number, ::Val{:Meeus})
    # Centuries to seconds conversion factor.
    cen2s = 1 / (36525.0 * 86400.0)

    # == Moon Position =====================================================================
    #
    # This is the same algorithm as in `moon_position_mod`. However, we must have access to
    # some internal variables. Hence, we copied the entire function for now.

    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000) / 36525

    # Moon's mean latitude referred to the mean equinox of date [deg].
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

    # Time derivatives of fundamental arguments [deg/cen].
    ∂L´_deg = @evalpoly(
        t_tdb,
        +481_267.881_234_21,
        2 * (-0.001_578_6),
        3 * (+1 / 538_841),
        4 * (-1 / 65_194_000)
    )

    ∂D_deg = @evalpoly(
        t_tdb,
        +445_267.111_403_4,
        2 * (-0.001_881_9),
        3 * (+1 / 545_868),
        4 * (-1 / 113_065_000)
    )

    ∂M_deg = @evalpoly(
        t_tdb,
        +35_999.050_290_9,
        2 * (-0.000_153_6),
        3 * (+1 / 24_490_000)
    )

    ∂M´_deg = @evalpoly(
        t_tdb,
        +477_198.867_505_5,
        2 * (+0.008_741_4),
        3 * (+1 / 69_699),
        4 * (-1 / 14_712_000)
    )

    ∂F_deg = @evalpoly(
        t_tdb,
        +483_202.017_523_3,
        2 * (-0.003_653_9),
        3 * (-1 / 3_526_000),
        4 * (+1 / 863_310_000)
    )

    ∂ϵ_deg = @evalpoly(t_tdb, -0.013_004_2, 2 * (-1.64e-7), 3 * (+5.04e-7))

    # Convert angles to [rad] and limit the angles between [0, 2π].
    L´ = mod2pi(deg2rad(L´))
    D  = mod2pi(deg2rad(D))
    M  = mod2pi(deg2rad(M))
    M´ = mod2pi(deg2rad(M´))
    F  = mod2pi(deg2rad(F))
    ϵ  = mod2pi(deg2rad(ϵ))
    A₁ = mod2pi(deg2rad(A₁))
    A₂ = mod2pi(deg2rad(A₂))
    A₃ = mod2pi(deg2rad(A₃))

    # Convert angular rates to [rad/cen].
    ∂L´ = deg2rad(∂L´_deg)
    ∂D  = deg2rad(∂D_deg)
    ∂M  = deg2rad(∂M_deg)
    ∂M´ = deg2rad(∂M´_deg)
    ∂F  = deg2rad(∂F_deg)

    # Term used to correct the arguments of the angle M that depends on the Earth's orbit
    # eccentricity around the Sun, and its time derivative.
    E  = @evalpoly(t_tdb, 1, -0.002_516, 0.000_007_4)
    E² = E * E
    ∂E = @evalpoly(t_tdb, -0.002_516, 2 * 0.000_007_4)

    # Compute the sum of the terms in the tables 47.A and 47.B and their time derivatives.

    tab = _TAB_47A
    num_terms = size(tab)[2]

    Σl  = 0.0
    Σr  = 0.0
    ∂Σl = 0.0
    ∂Σr = 0.0

    @inbounds for k in 1:num_terms
        aD  = tab[1, k]
        aM  = tab[2, k]
        aM´ = tab[3, k]
        aF  = tab[4, k]

        arg  = mod2pi(aD * D + aM * M + aM´ * M´ + aF * F)
        ∂arg = aD * ∂D + aM * ∂M + aM´ * ∂M´ + aF * ∂F

        E_corr, ∂E_corr = if ((aM == 1) || (aM == -1))
            E, ∂E
        elseif ((aM == 2) || (aM == -2))
            E², 2E * ∂E
        else
            one(E), zero(E)
        end

        sin_arg, cos_arg = sincos(arg)

        Σl += tab[5, k] * E_corr * sin_arg
        Σr += tab[6, k] * E_corr * cos_arg

        ∂Σl += tab[5, k] * (∂E_corr * sin_arg + E_corr * cos_arg * ∂arg)
        ∂Σr += tab[6, k] * (∂E_corr * cos_arg - E_corr * sin_arg * ∂arg)
    end

    tab = _TAB_47B
    num_terms = size(tab)[2]

    Σb  = 0.0
    ∂Σb = 0.0

    @inbounds for k in 1:num_terms
        aD  = tab[1, k]
        aM  = tab[2, k]
        aM´ = tab[3, k]
        aF  = tab[4, k]

        arg  = mod2pi(aD * D + aM * M + aM´ * M´ + aF * F)
        ∂arg = aD * ∂D + aM * ∂M + aM´ * ∂M´ + aF * ∂F

        E_corr, ∂E_corr = if ((aM == 1) || (aM == -1))
            E, ∂E
        elseif ((aM == 2) || (aM == -2))
            E², 2E * ∂E
        else
            one(E), zero(E)
        end

        sin_arg, cos_arg = sincos(arg)

        Σb  += tab[5, k] * E_corr * sin_arg
        ∂Σb += tab[5, k] * (∂E_corr * sin_arg + E_corr * cos_arg * ∂arg)
    end

    # Apply the corrections to the terms and their derivatives.
    ∂A₁ = deg2rad(131.849)
    ∂A₂ = deg2rad(479_264.290)
    ∂A₃ = deg2rad(481_266.484)

    sin_A₁, cos_A₁ = sincos(A₁)
    sin_A₂, cos_A₂ = sincos(A₂)
    sin_A₃, cos_A₃ = sincos(A₃)

    sin_L´_F, cos_L´_F     = sincos(L´ - F)
    sin_L´, cos_L´         = sincos(L´)
    sin_A₁_mF, cos_A₁_mF  = sincos(A₁ - F)
    sin_A₁_pF, cos_A₁_pF  = sincos(A₁ + F)
    sin_L´_mM´, cos_L´_mM´ = sincos(L´ - M´)
    sin_L´_pM´, cos_L´_pM´ = sincos(L´ + M´)

    Σl +=  3958sin_A₁ + 1962sin_L´_F + 318sin_A₂
    Σb += -2235sin_L´ + 382sin_A₃ + 175sin_A₁_mF + 175sin_A₁_pF +
           127sin_L´_mM´ - 115sin_L´_pM´

    ∂Σl += 3958cos_A₁ * ∂A₁ +
           1962cos_L´_F * (∂L´ - ∂F) +
           318cos_A₂ * ∂A₂

    ∂Σb += -2235cos_L´ * ∂L´ +
            382cos_A₃ * ∂A₃ +
            175cos_A₁_mF * (∂A₁ - ∂F) +
            175cos_A₁_pF * (∂A₁ + ∂F) +
            127cos_L´_mM´ * (∂L´ - ∂M´) -
            115cos_L´_pM´ * (∂L´ + ∂M´)

    # Convert to [rad].
    Σl = mod2pi(deg2rad(Σl / 1_000_000))
    Σb = mod2pi(deg2rad(Σb / 1_000_000))

    # Compute the Moon coordinates [rad] and [m].
    λ = mod2pi(L´ + Σl)
    β = Σb
    Δ = 385_000.56e3 + Σr

    sin_λ, cos_λ = sincos(λ)
    sin_β, cos_β = sincos(β)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    # == Moon Velocity =====================================================================

    # Convert rates to [rad/s] and [m/s].
    ∂λ = (∂L´ + deg2rad(∂Σl / 1_000_000)) * cen2s
    ∂β = deg2rad(∂Σb / 1_000_000) * cen2s
    ∂Δ = ∂Σr * cen2s
    ∂ϵ = deg2rad(∂ϵ_deg) * cen2s

    # Compute the Moon velocity vector in MOD [m/s].
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

    return v_moon_mod
end

function moon_velocity_mod(date_tdb::DateTime, ::Val{:Vallado})
    return moon_velocity_mod(datetime2julian(date_tdb), Val(:Vallado))
end

function moon_velocity_mod(jd_tdb::Number, ::Val{:Vallado})
    # Centuries to seconds conversion factor.
    cen2s = 1 / (36525.0 * 86400.0)

    # == Moon Position =====================================================================
    #
    # This is the same algorithm as in `moon_position_mod`. However, we must have access to
    # some internal variables. Hence, we copied the entire function for now.

    # Number of Julian centuries from J2000 epoch.
    t_tdb = (jd_tdb - JD_J2000) / 36525.0

    # Auxiliary computation to improve performance.
    sin1, cos1 = sincos(deg2rad(134.9 + 477_198.85t_tdb))
    sin2, cos2 = sincos(deg2rad(259.2 - 413_335.38t_tdb))
    sin3, cos3 = sincos(deg2rad(235.7 + 890_534.23t_tdb))
    sin4, cos4 = sincos(deg2rad(269.9 + 954_397.70t_tdb))
    sin5, cos5 = sincos(deg2rad(357.5 +  35_999.05t_tdb))
    sin6, cos6 = sincos(deg2rad(186.6 + 966_404.05t_tdb))

    # Ecliptic latitude of the Moon [deg].
    λₑ = 218.32 + 481_267.8813t_tdb + 6.29sin1 - 1.27sin2 + 0.66sin3 + 0.21sin4 - 0.19sin5 - 0.11sin6

    # Ecliptic longitude of the Moon [deg].
    sin_α₁, cos_α₁ = sincos(deg2rad( 93.3 + 483_202.03t_tdb))
    sin_α₂, cos_α₂ = sincos(deg2rad(228.2 + 960_400.87t_tdb))
    sin_α₃, cos_α₃ = sincos(deg2rad(318.3 +   6_003.18t_tdb))
    sin_α₄, cos_α₄ = sincos(deg2rad(217.6 - 407_332.20t_tdb))

    ϕₑ = 5.13sin_α₁ + 0.28sin_α₂ - 0.28sin_α₃ - 0.17sin_α₄

    # Parallax [deg].
    P = 0.9508 + 0.0518cos1 + 0.0095cos2 + 0.0078cos3 + 0.0028cos4

    # Obliquity of the ecliptic [deg].
    ϵ = @evalpoly(t_tdb, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)

    # Convert to radians and limit to the interval [0, 2π].
    λₑ = mod2pi(deg2rad(λₑ))
    ϕₑ = mod2pi(deg2rad(ϕₑ))
    P  = mod2pi(deg2rad(P))
    ϵ  = mod2pi(deg2rad(ϵ))

    # Compute the distance from Earth to the Moon [m].
    sin_P, cos_P = sincos(P)
    r = WGS84_ELLIPSOID.a / sin_P

    # Auxiliary variables.
    sin_λ, cos_λ = sincos(λₑ)
    sin_ϕ, cos_ϕ = sincos(ϕₑ)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    # == Moon Velocity =====================================================================

    # Time derivative of ecliptic latitude [deg/cen].
    ∂λₑ_deg = 481_267.8813 +
        6.29cos1  * deg2rad( 477_198.85) -
        1.27cos2  * deg2rad(-413_335.38) +
        0.66cos3  * deg2rad( 890_534.23) +
        0.21cos4  * deg2rad( 954_397.70) -
        0.19cos5  * deg2rad(  35_999.05) -
        0.11cos6  * deg2rad( 966_404.05)

    # Time derivative of ecliptic longitude [deg/cen].
    ∂ϕₑ_deg =
        5.13cos_α₁ * deg2rad( 483_202.03) +
        0.28cos_α₂ * deg2rad( 960_400.87) -
        0.28cos_α₃ * deg2rad(   6_003.18) -
        0.17cos_α₄ * deg2rad(-407_332.20)

    # Time derivative of parallax [deg/cen].
    ∂P_deg =
        -0.0518sin1 * deg2rad( 477_198.85) -
         0.0095sin2 * deg2rad(-413_335.38) -
         0.0078sin3 * deg2rad( 890_534.23) -
         0.0028sin4 * deg2rad( 954_397.70)

    # Time derivative of obliquity [deg/cen].
    ∂ϵ_deg = @evalpoly(t_tdb, -0.013_004_2, 2 * (-1.64e-7), 3 * (+5.04e-7))

    # Convert rates to [rad/s].
    ∂λ = deg2rad(∂λₑ_deg) * cen2s
    ∂ϕ = deg2rad(∂ϕₑ_deg) * cen2s
    ∂ϵ = deg2rad(∂ϵ_deg)  * cen2s
    ∂P = deg2rad(∂P_deg)  * cen2s

    # Distance rate [m/s]: r = a / sin(P), dr/dt = -a * cos(P) / sin²(P) * dP/dt.
    ∂r = -r * (cos_P / sin_P) * ∂P

    # Compute the Moon velocity vector represented in MOD (IAU-76/FK5 mean-equator,
    # mean-equinox of date) [m/s].
    v_moon_mod = SVector{3}(
        ∂r * cos_ϕ * cos_λ -
            r * ∂ϕ * sin_ϕ * cos_λ -
            r * ∂λ * cos_ϕ * sin_λ,
        ∂r * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ) -
            r * ∂ϵ * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ) -
            r * ∂ϕ * (cos_ϵ * sin_ϕ * sin_λ + sin_ϵ * cos_ϕ) +
            r * ∂λ * cos_ϵ * cos_ϕ * cos_λ,
        ∂r * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ) +
            r * ∂ϵ * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ) +
            r * ∂ϕ * (cos_ϵ * cos_ϕ - sin_ϵ * sin_ϕ * sin_λ) +
            r * ∂λ * sin_ϵ * cos_ϕ * cos_λ
    )

    return v_moon_mod
end
