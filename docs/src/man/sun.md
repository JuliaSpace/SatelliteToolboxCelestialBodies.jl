# Sun

```@meta
CurrentModule = SatelliteToolboxCelestialBodies
```

```@setup sun
using SatelliteToolboxCelestialBodies
```

The position [m] and velocity [m/s] of the Sun represented in the IAU-76/FK5 mean-equator,
mean-equinox of date (MOD) reference frame can be computed using the functions:

```julia
sun_position_mod(jd_tdb::Number) -> SVector{3, T}
sun_position_mod(date_tdb::Union{Date, DateTime}) -> SVector{3, Float64}
sun_velocity_mod(jd_tdb::Number) -> SVector{3, T}
sun_velocity_mod(date_tdb::Union{Date, DateTime}) -> SVector{3, Float64}
sun_state_mod(jd_tdb::Number) -> SVector{3, T}, SVector{3, T}
sun_state_mod(date_tdb::Union{Date, DateTime}) -> SVector{3, Float64}, SVector{3, Float64}
```

where the input epoch `jd_tdb` (Julian Day) or `date_tdb` must be represented in the
[Barycentric Dynamical Time (TDB)](https://en.wikipedia.org/wiki/Barycentric_Dynamical_Time).
The algorithm is described in **[1, pp. 277-279]**, and the velocity is obtained by
differentiating the position analytically.

The function `sun_state_mod` returns both the position and the velocity. Prefer it when
both quantities are required, since the computation is shared.

```@repl sun
sun_position_mod(DateTime(2006, 4, 2))

sun_velocity_mod(DateTime(2006, 4, 2))

sun_state_mod(date_to_jd(2006, 4, 2, 0, 0, 0))
```

!!! note

    The current civil time, as returned by `now()`, is neither TDB nor UTC. The difference
    between UTC and TDB is roughly 69 s, which is negligible for the accuracy of
    these models. Hence, `now(UTC)` can be used as the input epoch.

## References

- **[1]** Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*. 4th ed.
    Microcosm Press, Hawthorne, CA.
