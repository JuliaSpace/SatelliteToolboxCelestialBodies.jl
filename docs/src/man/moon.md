# Moon

```@meta
CurrentModule = SatelliteToolboxCelestialBodies
```

```@setup moon
using SatelliteToolboxCelestialBodies
```

The position [m] and velocity [m/s] of the Moon represented in the IAU-76/FK5 mean-equator,
mean-equinox of date (MOD) reference frame can be computed using the functions:

```julia
moon_position_mod(jd_tdb::Number[, model]) -> SVector{3, T}
moon_position_mod(date_tdb::Union{Date, DateTime}[, model]) -> SVector{3, Float64}
moon_velocity_mod(jd_tdb::Number[, model]) -> SVector{3, T}
moon_velocity_mod(date_tdb::Union{Date, DateTime}[, model]) -> SVector{3, Float64}
moon_state_mod(jd_tdb::Number[, model]) -> SVector{3, T}, SVector{3, T}
moon_state_mod(date_tdb::Union{Date, DateTime}[, model]) -> SVector{3, Float64}, SVector{3, Float64}
```

where the input epoch `jd_tdb` (Julian Day) or `date_tdb` must be represented in the
[Barycentric Dynamical Time (TDB)](https://en.wikipedia.org/wiki/Barycentric_Dynamical_Time).

The `model` selects the algorithm and must be `Val(:Meeus)` (default) or `Val(:Vallado)`:

- `Val(:Meeus)` uses the algorithm in **[2, ch. 47]** that provides an accuracy of 10 arcsec
  in the longitude and 4 arcsec in the latitude (the reference does not mention the
  timespan).
- `Val(:Vallado)` uses the algorithm in **[1, p. 288]** that is about 3 times faster than
  `Val(:Meeus)` but can lead to errors of 0.3° in longitude and 0.2° in latitude.

In both cases, the velocity is obtained by differentiating the position analytically. The
function `moon_state_mod` returns both the position and the velocity. Prefer it when both
quantities are required, since the computation is shared.

```@repl moon
moon_position_mod(DateTime(1994, 4, 28))

moon_position_mod(DateTime(1994, 4, 28), Val(:Vallado))

moon_velocity_mod(DateTime(1994, 4, 28))

moon_state_mod(date_to_jd(1994, 4, 28, 0, 0, 0), Val(:Vallado))
```

!!! note

    The current civil time, as returned by `now()`, is neither TDB nor UTC. The difference
    between UTC and TDB is roughly 69 s, which is negligible for the accuracy of
    these models. Hence, `now(UTC)` can be used as the input epoch.

## References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
- **[2]** Meeus, J. (1998). *Astronomical Algorithms*. 2nd ed. Willmann-Bell, Inc,
    Richmond, VA.
