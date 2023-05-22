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
function sun_position_mod(jd_tdb::Number)
function sun_position_mod(date_tdb::DateTime)
```

where the input time `jd_tdb` (Julian Day) or `date_tdb` must be represented in the
[Barycentric Dynamical Time (TDB)](https://en.wikipedia.org/wiki/Barycentric_Dynamical_Time).

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
function sun_velocity_mod(jd_tdb::Number)
function sun_velocity_mod(date_tdb::DateTime)
```

where the input time `jd_tdb` (Julian Day) or `date_tdb` must be represented in the
[Barycentric Dynamical Time (TDB)](https://en.wikipedia.org/wiki/Barycentric_Dynamical_Time).

> **Note**
> This algorithm was obtained by differentiating the Sun position equations in **[1]**.

```julia
julia> sun_velocity_mod(now())
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -25645.525387897742
  13231.300593568181
   5735.626095163374
```

### Moon

We can compute the Moon position represented in the Mean-Of-Date (MOD) reference frame
**[1, 2]** using the functions:

``` julia
function moon_position_mod(jd_tdb::Number[, model]) -> SVector{3, Float64}
function moon_position_mod(date_tdb::DateTime[, model]) -> SVector{3, Float64}
```

where the input time `jd_tdb` (Julian Day) or `date_tdb` must be represented in the
Barycentric Dynamical Time (TDB).

The `model` must be `Val(:Meeus)` or `Val(:Vallado)`. `Val(:Meeus)` uses the algorithm in
**[2, p. 337]** that provides an accuracy of 10" in the longitude and 4" in the latitude
(the reference does not mention the timespan). `Val(:Vallado)` uses the algorithm in
**[1, p. 288]** that is 10x faster than `Val(:Meeus)` but can lead to errors of 0.3° in
longitude and 0.2° in latitude.

```julia
julia> moon_position_mod(now())
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -4.992612797700085e7
  3.48593091076279e8
  1.864034978650991e8

julia> moon_position_mod(now(), Val(:Vallado))
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -4.991011671868989e7
  3.481318482554912e8
  1.8647115876567587e8
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
- **[2]** **Meeus, J** (1998). *Astronomical algorithms*. **Willmann-Bell, Inc**, Richmond, VA.
