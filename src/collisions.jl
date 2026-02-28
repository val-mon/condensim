function compute_next_pos!(mol, dt)
    acc = zeros(3)
    mol.velocity = mol.velocity + acc * dt
    mol.pos = mol.pos + mol.velocity * dt
end

function are_colliding(mol1::Molecule, mol2::Molecule)::Bool
    return norm(mol1.pos .- mol2.pos) <= mol1.radius + mol2.radius
end

function resolve_collision!(mol1::Molecule, mol2::Molecule)
    # compute unit vector of the center from m1 to m2
    diff = mol2.pos .- mol1.pos
    dist = norm(diff)
    n̂ = diff ./ dist

    # compute scalar product of relative velocity along n̂
    Δv = dot(mol1.velocity .- mol2.velocity, n̂)

    # if molecules are already moving away – no shock to resolve
    if Δv <= 0
        return
    end

    # otherwise compute new velocities
    mol1.velocity = mol1.velocity .- (2 * mol2.mass / (mol1.mass + mol2.mass)) .* Δv .* n̂
    mol2.velocity = mol2.velocity .+ (2 * mol1.mass / (mol1.mass + mol2.mass)) .* Δv .* n̂
end

function reflect_walls!(mol::Molecule, box)
    for i in 1:3
        if mol.pos[i] - mol.radius < -box
            mol.pos[i] = -box + mol.radius
            mol.velocity[i] = abs(mol.velocity[i])
        elseif mol.pos[i] + mol.radius > box
            mol.pos[i] = box - mol.radius
            mol.velocity[i] = -abs(mol.velocity[i])
        end
    end
end
