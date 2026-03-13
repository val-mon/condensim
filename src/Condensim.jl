module Condensim

using Plots
using Random
using LinearAlgebra
using Statistics

include("types.jl")
include("physics.jl")
include("collisions.jl")
include("simulation.jl")

export Domain, volume, reflect_walls!
export Molecule, resolve_collision!, are_colliding

function main()
    simulation()
end

end
