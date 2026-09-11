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

############################################################################################
#                                        Constants                                         #
############################################################################################

# Conversion factor from Julian centuries to seconds [s/century].
const _CENTURY_TO_SECONDS = 1 / (36525 * 86400)

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./moon/constants.jl")
include("./moon/moon.jl")

include("./sun.jl")

end # module SatelliteToolboxCelestialBodies
