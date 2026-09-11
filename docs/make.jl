using Documenter
using SatelliteToolboxCelestialBodies

makedocs(;
    modules = [SatelliteToolboxCelestialBodies],
    format = Documenter.HTML(;
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteToolboxCelestialBodies.jl/stable/",
    ),
    sitename = "SatelliteToolboxCelestialBodies.jl",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Celestial Bodies" => [
            "Sun"  => "man/sun.md",
            "Moon" => "man/moon.md",
        ],
        "Library" => "lib/library.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl.git", target = "build"
)
