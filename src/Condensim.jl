module Condensim

using Plots
using Random

include("types.jl")
include("physics.jl")
include("mouvement.jl")

function main()
    simulate_mouvement()
end

end
