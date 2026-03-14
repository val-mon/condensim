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

function resolve_collisions!(molecules)
    for i in 1:length(molecules)
        for j in i+1:length(molecules)
            if are_colliding(molecules[i], molecules[j])
                resolve_collision!(molecules[i], molecules[j])
            end
        end
    end
end

function reflect_walls!(mol::Molecule, domain::Domain)
    half = [domain.Lx, domain.Ly, domain.Lz] ./ 2
    for i in 1:3
        lo = -half[i] + mol.radius
        hi =  half[i] - mol.radius
        if mol.pos[i] < lo
            mol.pos[i] = 2*lo - mol.pos[i]
            mol.velocity[i] = abs(mol.velocity[i])
        elseif mol.pos[i] > hi
            mol.pos[i] = 2*hi - mol.pos[i]
            mol.velocity[i] = -abs(mol.velocity[i])
        end
    end
end
