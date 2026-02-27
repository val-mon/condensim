mutable struct Vec3
    x::Float64
    y::Float64
    z::Float64
end

function Base.show(io::IO, v::Vec3)
    print(io, "($(v.x), $(v.y), $(v.z))")
end


mutable struct Molecule
    # mass [kg]
    mass::Float64

    # sphere radius [m]
    radius::Float64

    # spatial position
    pos::Vec3
    
    # velocity [m/s]
    velocity::Vec3

    # chemical formula
    formula::String
end

function Base.show(io::IO, m::Molecule)
    println(io, "molecule : $(m.formula)")
    println(io, "  mass     [kg]  : $(m.mass)")
    println(io, "  radius   [m]   : $(m.radius)")
    println(io, "  position [m]   : $(m.pos)")
    println(io, "  velocity [m/s] : $(m.velocity)")
end
