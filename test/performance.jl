## Description #############################################################################
#
# Tests related to code quality, type inference, allocations, and automatic differentiation.
#
############################################################################################

# == Code Quality ==========================================================================

@testset "Aqua.jl" begin
    Aqua.test_all(SatelliteToolboxCelestialBodies)
end

if VERSION >= v"1.12"
    @warn "JET.jl test skipped on Julia 1.12+ due to MethodTableView incompatibility"
else
    @testset "JET.jl" begin
        JET.test_package(
            SatelliteToolboxCelestialBodies;
            toplevel_logger = nothing,
            target_modules = (SatelliteToolboxCelestialBodies,),
        )
    end
end

# == Type Inference and Allocations ========================================================

# Measure the allocations of `f(args...)` after a warm-up call.
function _allocations(f, args...)
    f(args...)
    return @allocated f(args...)
end

@testset "Type Inference and Allocations" begin
    jd_tdb   = date_to_jd(2000, 6, 15, 12, 0, 0)
    date_tdb = jd_tdb |> julian2datetime
    day_tdb  = Date(2000, 6, 15)

    for f in (sun_position_mod, sun_velocity_mod, sun_state_mod)
        @test @inferred(f(jd_tdb)) == f(jd_tdb)
        @test @inferred(f(date_tdb)) == f(date_tdb)
        @test @inferred(f(day_tdb)) == f(day_tdb)
        @test _allocations(f, jd_tdb) == 0
    end

    for f in (moon_position_mod, moon_velocity_mod, moon_state_mod)
        @test @inferred(f(jd_tdb)) == f(jd_tdb)
        @test @inferred(f(date_tdb)) == f(date_tdb)
        @test @inferred(f(day_tdb)) == f(day_tdb)

        for model in (Val(:Meeus), Val(:Vallado))
            @test @inferred(f(jd_tdb, model)) == f(jd_tdb, model)
            @test @inferred(f(date_tdb, model)) == f(date_tdb, model)
            @test _allocations(f, jd_tdb, model) == 0
        end
    end
end

# == Element Type ==========================================================================

@testset "Element Type Promotion" begin
    jd_tdb = date_to_jd(2000, 6, 15, 12, 0, 0)

    # `Float32` inputs must be promoted to `Float64`.
    for f in (sun_position_mod, sun_velocity_mod)
        @test f(Float32(jd_tdb)) isa SVector{3, Float64}
        @test f(round(Int, jd_tdb)) isa SVector{3, Float64}
        @test f(big(jd_tdb)) isa SVector{3, BigFloat}
    end

    for f in (moon_position_mod, moon_velocity_mod), model in (Val(:Meeus), Val(:Vallado))
        @test f(Float32(jd_tdb), model) isa SVector{3, Float64}
        @test f(round(Int, jd_tdb), model) isa SVector{3, Float64}
        @test f(big(jd_tdb), model) isa SVector{3, BigFloat}
    end
end

# == Automatic Differentiation =============================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# The analytical velocities must match the derivative of the position functions obtained
# with ForwardDiff.jl. Notice that the derivative w.r.t. the Julian day is in [m/day].
#
############################################################################################

@testset "Automatic Differentiation" begin
    jd_start = date_to_jd(1950, 1, 1, 0, 0, 0)
    jd_stop  = date_to_jd(2100, 1, 1, 0, 0, 0)

    for _ in 1:100
        jd_tdb = jd_start + rand() * (jd_stop - jd_start)

        v_ad = ForwardDiff.derivative(sun_position_mod, jd_tdb) / 86400
        @test v_ad ≈ sun_velocity_mod(jd_tdb) rtol = 1e-10

        for model in (Val(:Meeus), Val(:Vallado))
            v_ad =
                ForwardDiff.derivative(jd -> moon_position_mod(jd, model), jd_tdb) / 86400
            @test v_ad ≈ moon_velocity_mod(jd_tdb, model) rtol = 1e-10
        end
    end
end
