using Test

using LinearAlgebra
using ReferenceFrameRotations
using SatelliteToolboxCelestialBodies

@testset "Sun" verbose = true begin
    include("./sun.jl")
end

@testset "Moon" verbose = true begin
    include("./moon.jl")
end
