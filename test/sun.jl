# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the position of the Sun.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Blanco, Manuel Jesus, Milidonis, Kypros, Bonanos, Aristides. Updating the PSA sun
#       position algorithm. Solar Energy, vol.212, Elsevier BV,2020-12.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/sun.jl
# ==========================================================================================

# Function sun_position_mod
# ------------------------------------------------------------------------------------------

############################################################################################
#                                       Test Results
############################################################################################
#
# Example 5-1: Finding the Sun position vector [1, p. 280]
#
#   Using UTC = April 2, 2006, 00:00, and jd_utc ≈ jd_tdb, one gets:
#
#           | 0.9771945 |
#   r_mod = | 0.1924424 | AU
#           | 0.0834308 |
#
############################################################################################

@testset "Sun Position" begin
    s_mod = sun_position_mod(DateTime("2006-04-02T00:00:00")) / ASTRONOMICAL_UNIT

    @test s_mod[1] ≈ 0.9771945 atol = 2e-6
    @test s_mod[2] ≈ 0.1924424 atol = 2e-6
    @test s_mod[3] ≈ 0.0834308 atol = 2e-6
end

# Function sun_velocity_mode
# ------------------------------------------------------------------------------------------

############################################################################################
#                                       Test Results
############################################################################################
#
# The time-derivative of the Sun position vector computed by the function `sun_velocity_mod`
# is compared to the numerical differentiation of the Sun position vector using the function
# `sun_position_mod`, which is already validated.
#
############################################################################################

@testset "Function sun_velocity_i" begin
    jd_start = date_to_jd(1950, 1, 1, 0, 0, 0)
    jd_stop  = date_to_jd(2019, 1, 1, 0, 0, 0)

    for _ in 1:100
        jd_tdb = rand(jd_start:jd_stop)
        Δt     = 0.1
        s_t₁   = sun_position_mod(jd_tdb |> julian2datetime)
        s_t₂   = sun_position_mod(jd_tdb + Δt / 86400 |> julian2datetime)
        v_n    = (s_t₂ - s_t₁) / Δt
        v      = sun_velocity_mod(jd_tdb |> julian2datetime)

        @test norm(v - v_n) / norm(v) * 100 < 0.055
    end
end

############################################################################################
#                                       Test Results
############################################################################################
#
# Finding the local sun position in local frame using the following data:
#   JD = 2.458985282701389e6
#   Latitude = 37.1 °N
#   Longitude = -2.36 °W
#
# Calling the function as:
#   sun_position_el(2.458985282701389e6, 37.1, -2.36, 'a')
#
# Must result in :
#   RA ~ 53.04443267321162 °
#   DEC ~ 19.108142467612954 °
#   HA ~ 100.31967018417757 °
#   ZEN ~ 86.42495381753338 °
#   AZM ~ 291.337579668121 °
#   SUN ~ [-0.9296, 0.3632, 0.0624]
#
# Accuracy: 50 arcsecs ~ 0.0138889 ° [2, p. 341]
#
############################################################################################

@testset "Sun Position Local" begin
    JD = 2.458985282701389e6
    Latitude = 37.1
    Longitude = -2.36

    # Tolerances
    tol_ang = 1.38889e-2
    tol_vec = 1e-2

    (ra, dec, ha, zen, az, sun) = sun_position_el(JD, Latitude, Longitude, 'a')

    @test ra ≈ 53.04443267321162 atol = tol_ang
    @test dec ≈ 19.108142467612954 atol = tol_ang
    @test ha ≈ 100.31967018417757 atol = tol_ang
    @test zen ≈ 86.42495381753338 atol = tol_ang
    @test az ≈ 291.337579668121 atol = tol_ang
    @test sun ≈ [-0.9296, 0.3632, 0.0624] atol = tol_vec
end
