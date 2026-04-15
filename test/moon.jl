## Description #############################################################################
#
# Tests related to the position and velocity of the Moon.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorn, CA, USA.
#
# [2] Meeus, J (1998). Astronomical algorithms. Willmann-Bell, Inc, Richmond, VA.
#
############################################################################################

# == File: ./src/moon/moon.jl ==============================================================

# -- Function moon_position_mod ------------------------------------------------------------

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# Example 47.a: Calculate the geocentric longitude, latitude, and equatorial horizontal
#               parallax of the Moon for 1992 April 12, at 0h TD.
#
#   In the ecliptic plane w.r.t. the mean equator of date:
#
#       λ = 133°.162655 (longitude)
#       β = -3°.229126 (latitude)
#       Δ = 368409.7 km
#
############################################################################################

@testset "Moon Position (Meeus)" begin
    date_tdb = DateTime("1992-04-12T00:00:00")

    r_moon_mod = moon_position_mod(date_tdb)

    # We need to compute the mean ecliptic to obtain the variables in the example.
    t_tdb = (datetime2julian(date_tdb) - JD_J2000) / 36525
    ϵ = @evalpoly(t_tdb, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)
    ϵ = mod2pi(deg2rad(ϵ))

    r_moon_ecl = angle_to_dcm(ϵ, :X) * r_moon_mod

    λ = atand(r_moon_ecl[2], r_moon_ecl[1])
    β = atand(r_moon_ecl[3], √(r_moon_ecl[1]^2 + r_moon_ecl[2]^2))
    Δ = norm(r_moon_ecl) / 1000

    @test λ ≈ 133.162655 atol = 1e-6
    @test β ≈ -3.229126  atol = 1e-6
    @test Δ ≈ 368409.7   atol = 1e-1

    # Test overloads.
    r_moon2_mod = moon_position_mod(date_tdb |> datetime2julian)

    @test r_moon2_mod == r_moon_mod
end

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# Example 5-3: Finding the Moon position vector [1, p. 289]
#
#   In April 28, 1994, 00:00 TDB, one gets:
#
#           | -134_240.626 |
#   r_mod = | -311_571.590 | km
#           | -126_693.785 |
#
############################################################################################

@testset "Moon Position (Vallado)" begin
    date_tdb = DateTime("1994-04-28T00:00:00")
    r_mod = moon_position_mod(date_tdb, Val(:Vallado))

    @test r_mod[1] / 1000 ≈ -134_240.626 atol = 1e-3
    @test r_mod[2] / 1000 ≈ -311_571.590 atol = 1e-3
    @test r_mod[3] / 1000 ≈ -126_693.785 atol = 1e-3

    # Test overloads.
    r_mod2 = moon_position_mod(date_tdb |> datetime2julian, Val(:Vallado))

    @test r_mod2 == r_mod
end

@testset "Moon Position (Overloads)" begin
    jd_tdb   = date_to_jd(2000, 6, 15, 12, 0, 0)
    date_tdb = julian2datetime(jd_tdb)

    # Default model is Meeus.
    @test moon_position_mod(jd_tdb) == moon_position_mod(jd_tdb, Val(:Meeus))
    @test moon_position_mod(date_tdb) == moon_position_mod(date_tdb, Val(:Meeus))

    # DateTime and Julian Day must give the same result.
    @test moon_position_mod(date_tdb, Val(:Meeus)) == moon_position_mod(jd_tdb, Val(:Meeus))
    @test moon_position_mod(date_tdb, Val(:Vallado)) == moon_position_mod(jd_tdb, Val(:Vallado))
end

# -- Function moon_velocity_mod -----------------------------------------------------------

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# The time-derivative of the Moon position vector computed by the function
# `moon_velocity_mod` is compared to the numerical differentiation of the Moon position
# vector using the function `moon_position_mod`, which is already validated.
#
############################################################################################

@testset "Moon Velocity (Meeus)" begin
    jd_start = date_to_jd(1950, 1, 1, 0, 0, 0)
    jd_stop  = date_to_jd(2019, 1, 1, 0, 0, 0)

    for _ in 1:100
        jd_tdb = rand(jd_start:jd_stop)
        Δt     = 0.1
        r_t₁   = moon_position_mod(jd_tdb |> julian2datetime)
        r_t₂   = moon_position_mod(jd_tdb + Δt / 86400 |> julian2datetime)
        v_n    = (r_t₂ - r_t₁) / Δt
        v      = moon_velocity_mod(jd_tdb |> julian2datetime)

        @test norm(v - v_n) / norm(v) * 100 < 0.055
    end
end

@testset "Moon Velocity (Vallado)" begin
    jd_start = date_to_jd(1950, 1, 1, 0, 0, 0)
    jd_stop  = date_to_jd(2019, 1, 1, 0, 0, 0)

    for _ in 1:100
        jd_tdb = rand(jd_start:jd_stop)
        Δt     = 0.1
        r_t₁   = moon_position_mod(jd_tdb |> julian2datetime, Val(:Vallado))
        r_t₂   = moon_position_mod(jd_tdb + Δt / 86400 |> julian2datetime, Val(:Vallado))
        v_n    = (r_t₂ - r_t₁) / Δt
        v      = moon_velocity_mod(jd_tdb |> julian2datetime, Val(:Vallado))

        @test norm(v - v_n) / norm(v) * 100 < 0.055
    end
end

@testset "Moon Velocity (Overloads)" begin
    jd_tdb   = date_to_jd(2000, 6, 15, 12, 0, 0)
    date_tdb = julian2datetime(jd_tdb)

    # Default model is Meeus.
    @test moon_velocity_mod(jd_tdb) == moon_velocity_mod(jd_tdb, Val(:Meeus))
    @test moon_velocity_mod(date_tdb) == moon_velocity_mod(date_tdb, Val(:Meeus))

    # DateTime and Julian Day must give the same result.
    @test moon_velocity_mod(date_tdb, Val(:Meeus)) == moon_velocity_mod(jd_tdb, Val(:Meeus))
    @test moon_velocity_mod(date_tdb, Val(:Vallado)) == moon_velocity_mod(jd_tdb, Val(:Vallado))
end
