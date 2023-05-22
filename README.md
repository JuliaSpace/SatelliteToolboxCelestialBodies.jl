SatelliteToolboxCelestialBodies.jl
==================================

This package contains functions to compute the position and velocity of some celestial
bodies for the **SatelliteToolbox.jl** ecosystem.

## Installation

```julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxCelestialBodies.jl")
```

## Usage

### Sun

We can compute the Sun position represented in the Mean-Of-Date (MOD) reference frame
**[1]** using the functions:

```julia
function sun_position_mod(jd::Number)
function sun_position_mod(date::DateTime)
```

The input date must be represented in [UT1](https://en.wikipedia.org/wiki/Universal_Time).

```julia
julia> sun_position_mod(now())
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 7.281649894711235e10
 1.2182511371727788e11
 5.2809968734836815e10
```

We can also compute the Sun velocity represented in MOD frame as measured by an observer in
the same frame:

```julia
function sun_velocity_mod(jd::Number)
function sun_velocity_mod(date::DateTime)
```

The input date must be represented in [UT1](https://en.wikipedia.org/wiki/Universal_Time).

> **Note**
> This algorithm was obtained by differentiating the Sun position equations in **[1]**.

```julia
julia> sun_velocity_mod(now())
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -25645.525387897742
  13231.300593568181
   5735.626095163374
```

## Rationale

The packages in [JuliaAstro](https://github.com/JuliaAstro) provide the same functionality
with usually more precision than the algorithms here. However, our goal is to build an
attitude and orbit control subsystem written in Julia. In this case, the footprint of each
package matters. Hence, we created this small package to contain the necessary functions
since the extensive feature list in the other packages is unnecessary here.

## References

- **[1]** **Vallado, D. A** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorn, CA, USA.
