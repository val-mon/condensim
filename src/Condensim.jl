module Condensim

using Plots
using Random
using LinearAlgebra

include("types.jl")
include("physics.jl")
include("collisions.jl")
include("simulation.jl")

function main()
    simulation()
end

end
