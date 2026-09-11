# SatelliteToolboxCelestialBodies.jl

This package contains functions to compute the position and velocity of celestial bodies
for the **SatelliteToolbox.jl** ecosystem. Currently, the following bodies are available:

- Sun, using the algorithm in **[1]**; and
- Moon, using the algorithms in **[1]** and **[2]**.

## Installation

This package can be installed using:

```julia
julia> using Pkg
julia> Pkg.add("SatelliteToolboxCelestialBodies")
```

## Rationale

The packages in [JuliaAstro](https://github.com/JuliaAstro) provide the same functionality
with usually more precision than the algorithms here. However, our goal is to build an
attitude and orbit control subsystem written in Julia. In this case, the footprint of each
package matters. Hence, we created this small package to contain the necessary functions
since the extensive feature list in the other packages is unnecessary here.

## References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
- **[2]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.
