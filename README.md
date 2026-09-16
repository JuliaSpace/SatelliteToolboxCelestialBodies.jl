<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxCelestialBodies.jl"><br>
  <small><i>This package is part of the <a href="https://github.com/JuliaSpace/SatelliteToolbox.jl">SatelliteToolbox.jl</a> ecosystem.</i></small>
</p>

# SatelliteToolboxCelestialBodies.jl

[![CI](https://img.shields.io/github/actions/workflow/status/JuliaSpace/SatelliteToolboxCelestialBodies.jl/ci.yml?style=flat-square&logo=githubactions&logoColor=white&labelColor=475569&label=CI)](https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl/actions/workflows/ci.yml)
[![Codecov](https://img.shields.io/codecov/c/github/JuliaSpace/SatelliteToolboxCelestialBodies.jl?token=CONQMSI4JD&style=flat-square&logo=codecov&logoColor=white&labelColor=475569)](https://codecov.io/gh/JuliaSpace/SatelliteToolboxCelestialBodies.jl)
[![docs-stable](https://img.shields.io/badge/docs-stable-16A34A?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-stable-url]
[![docs-dev](https://img.shields.io/badge/docs-dev-D97706?style=flat-square&logo=gitbook&logoColor=white&labelColor=475569)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495D1?style=flat-square&labelColor=475569)](https://github.com/invenia/BlueStyle)
[![License](https://img.shields.io/github/license/JuliaSpace/SatelliteToolboxCelestialBodies.jl?style=flat-square&logo=readme&logoColor=white&labelColor=475569&color=0284C7)](https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl/blob/main/LICENSE.txt)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.11262848-DB2777?style=flat-square&logo=doi&logoColor=white&labelColor=475569)](https://zenodo.org/doi/10.5281/zenodo.11262848)

This package contains functions to compute the position and velocity of some celestial
bodies for the **SatelliteToolbox.jl** ecosystem.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("SatelliteToolboxCelestialBodies")
```

## Usage

The position [m] and velocity [m/s] of the Sun and Moon represented in the IAU-76/FK5
mean-equator, mean-equinox of date (MOD) reference frame can be computed using the
functions `sun_position_mod`, `sun_velocity_mod`, `sun_state_mod`, `moon_position_mod`,
`moon_velocity_mod`, and `moon_state_mod`. The input epoch is a Julian Day, a `Date`, or a
`DateTime` in the
[Barycentric Dynamical Time (TDB)](https://en.wikipedia.org/wiki/Barycentric_Dynamical_Time).

```julia
julia> sun_position_mod(now(UTC))
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 7.281649894711235e10
 1.2182511371727788e11
 5.2809968734836815e10

julia> r_moon_mod, v_moon_mod = moon_state_mod(now(UTC), Val(:Vallado));
```

> [!NOTE]
> The difference between UTC and TDB is roughly 69 s, which is negligible for the accuracy
> of these models.

See the [package documentation][docs-stable-url] for more details.

## Rationale

The packages in [JuliaAstro](https://github.com/JuliaAstro) provide the same functionality
with usually more precision than the algorithms here. However, our goal is to build an
attitude and orbit control subsystem written in Julia. In this case, the footprint of each
package matters. Hence, we created this small package to contain the necessary functions
since the extensive feature list in the other packages is unnecessary here.

## References

- **[1]** **Vallado, D. A.** (2013). *Fundamentals of Astrodynamics and Applications*. 4th
  ed. **Microcosm Press**, Hawthorne, CA.
- **[2]** **Meeus, J.** (1998). *Astronomical Algorithms*. 2nd ed. **Willmann-Bell, Inc**,
  Richmond, VA.

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolboxCelestialBodies.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolboxCelestialBodies.jl/stable
