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
#   JD = 2460065.945463
#   Latitude = 51.0 °N
#   Longitude = 10.0 °E
#   Pressure = 1.0 atm
#   Temperature = 20.0 °C
#   Algorithm = 5
#
# Calling the function as:
#   sun_position_el(2460065.945463, 51.0, 10.0, 1.0, 20.0, 'e', '5')
#
# Must result in :
#   RA = 02h 34m 09s; ~ 38.5375 °
#   DEC = +15° 06’ 42”; ~ 15.111667 °
#
############################################################################################

@testset "Sun Position Local" begin
    JD = 2460065.945463
    Latitude = 51.0
    Longitude = 10.0
    Pressure = 1.0
    Temperature = 20.0

    # Tolerances for various algorithms
    alg = ('1', '2', '3', '4', '5')
    tol = (2e0, 2e0, 1e0, 1e0, 1e0)

    for algorithm in 1:5
        s_eq = sun_position_el(JD, Latitude, Longitude, Pressure, Temperature, 'e', alg[algorithm])

        # Convert outputs to degrees
        ra = s_eq[1]*180/π
        dec = s_eq[2]*180/π

        @test ra ≈ 38.5375 atol = tol[algorithm]
        @test dec ≈ 15.111667 atol = tol[algorithm]
    end
end
