module Condensim

# simulation libs
using Plots

include("types.jl")
include("physics.jl")

function main()
    # helium
    He = Molecule(
        4.0026 * PHYSICS.u,
        31e-12,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        "He"
    )
    println(He)

    # neon
    Ne = Molecule(
        20.1797 * PHYSICS.u,
        38e-12,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        "Ne"
    )
    println(Ne)

    # nitrogen
    N2 = Molecule(
        28.0134 * PHYSICS.u,
        65e-12,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        "N2"
    )
    println(N2)

    # oxygen
    O2 = Molecule(
        31.9988 * PHYSICS.u,
        60e-12,
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.0, 0.0, 0.0),
        "O2"
    )
    println(O2)
end

end
