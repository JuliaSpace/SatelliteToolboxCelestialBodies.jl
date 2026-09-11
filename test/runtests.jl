## Description #############################################################################
#
# Run the test suite of SatelliteToolboxCelestialBodies.jl.
#
############################################################################################

using Test

using LinearAlgebra
using ReferenceFrameRotations
using SatelliteToolboxCelestialBodies
using StaticArrays

@testset "Sun" verbose = true begin
    include("./sun.jl")
end

@testset "Moon" verbose = true begin
    include("./moon.jl")
end

if isempty(VERSION.prerelease)
    using Aqua
    using ForwardDiff
    using JET

    @testset "Performance" verbose = true begin
        include("./performance.jl")
    end
else
    @warn "Performance checks not guaranteed to work on julia-nightly, skipping tests"
end
