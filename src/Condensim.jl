module Condensim

using Plots
using Random
using LinearAlgebra
using Statistics

include("types.jl")
include("physics.jl")
include("utils.jl")
include("simulation.jl")
include("shannon.jl")
include("thermal_wall.jl")
include("lennard_jones.jl")

export Domain, volume, reflect_walls!
export Molecule, resolve_collision!, are_colliding

function main()
    launch_simulation()
end

end
