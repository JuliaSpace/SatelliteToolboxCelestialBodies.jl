SatelliteToolboxCelestialBodies.jl Changelog
============================================

Version 1.1.0
-------------

- ![Feature][badge-feature] The functions `moon_velocity_mod` (PR [#8][gh-pr-8]),
  `sun_state_mod`, and `moon_state_mod` were added. The latter two return the position and
  velocity at once, sharing the computation.
- ![Feature][badge-feature] All functions now accept a `Date` as the input epoch.
- ![Feature][badge-feature] The Moon functions now throw an `ArgumentError` with a helpful
  message if the selected model is not supported.
- ![Enhancement][badge-enhancement] The algorithms are now generic in the numeric type of
  the input, propagating types such as `BigFloat` and automatic differentiation numbers.
  `Float32` inputs are promoted to `Float64`.
- ![Enhancement][badge-enhancement] The sums over the Meeus tables are now generated at
  compile time with literal multipliers, and the sines and cosines of the arguments are
  composed from those of the fundamental arguments using the angle addition formulas
  instead of one `sincos` call per term. The Meeus model is now roughly 4 times faster.
- ![Enhancement][badge-enhancement] The duplicated code between the position and velocity
  functions was removed.
- ![Enhancement][badge-enhancement] The package now has a documentation site and a
  precompilation workload.
- ![Enhancement][badge-enhancement] The docstrings were rewritten to document units,
  reference frames, time scales, and the element type of the results.
- ![Bugfix][badge-bugfix] The time derivative of the Sun distance used `cos(2Ms)` instead of
  `sin(2Ms)`, leading to errors of up to 8 m/s in the Sun velocity.
- ![Bugfix][badge-bugfix] The sign of the T² term in the eccentricity correction factor `E`
  of the Meeus model was wrong. The error is negligible near J2000 but grows quadratically
  with time.
- ![Info][badge-info] The test suite now checks the code quality with Aqua.jl and JET.jl,
  the type inference and allocations of all public methods, and compares the analytical
  velocities with derivatives obtained with ForwardDiff.jl.
- ![Info][badge-info] The package is compatible with SatelliteToolboxBase.jl 1 and 2.

Version 1.0.1
-------------

- ![Enhancement][badge-enhancement] The package now supports Zygote 0.7. (PR [#4][gh-pr-4])

Version 1.0.0
-------------

- ![Info][badge-info] We dropped support for Julia 1.6. This version only supports the
  current Julia version and v1.10 (LTS).
- ![Info][badge-info] This version does not have breaking changes. We bump the version to
  1.0.0 because we now consider the API stable.

Version 0.1.2
-------------

- ![Enhancement][badge-enhancement] Documentation update and minor source-code improvements.

Version 0.1.1
-------------

- ![Enhancement][badge-enhancement] We updated the dependency compatibility bounds.

Version 0.1.0
-------------

- Initial version.
  - This version was based on the code in **SatelliteToolbox.jl**.

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

[gh-pr-4]: https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl/pull/4
[gh-pr-8]: https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl/pull/8
