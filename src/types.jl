mutable struct Molecule
    # mass [kg]
    mass::Float64

    # sphere radius [m]
    radius::Float64

    # spatial position
    pos::Vector{Float64}
    
    # velocity [m/s]
    velocity::Vector{Float64}

    # chemical formula
    formula::String
end

function Base.copy(m::Molecule)
    Molecule(m.mass, m.radius, copy(m.pos), copy(m.velocity), m.formula)
end

function Base.show(io::IO, m::Molecule)
    println(io, "molecule : $(m.formula)")
    println(io, "  mass     [kg]  : $(m.mass)")
    println(io, "  radius   [m]   : $(m.radius)")
    println(io, "  position [m]   : $(m.pos)")
    println(io, "  velocity [m/s] : $(m.velocity)")
end
