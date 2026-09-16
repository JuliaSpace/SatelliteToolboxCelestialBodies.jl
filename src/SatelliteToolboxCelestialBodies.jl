"""
    module SatelliteToolboxCelestialBodies

Compute the position and velocity of celestial bodies (Sun and Moon) for the
**SatelliteToolbox.jl** ecosystem.
"""
module SatelliteToolboxCelestialBodies

using Reexport
using StaticArrays

@reexport using Dates
@reexport using SatelliteToolboxBase

import PrecompileTools

############################################################################################
#                                        Constants                                         #
############################################################################################

# Conversion factor from rates per Julian century to rates per second [century/s].
const _CENTURIES_PER_SECOND = 1 / (36525 * 86400)

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./moon/constants.jl")
include("./moon/moon.jl")

include("./math.jl")
include("./sun.jl")

include("./precompile.jl")

end # module SatelliteToolboxCelestialBodies
