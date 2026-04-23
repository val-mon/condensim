const MOL_TYPES = [
    Molecule(4.0026 * PHYSICS.u, 1.1e-10, zeros(3), zeros(3), "He"),
    Molecule(20.1797 * PHYSICS.u, 38e-12, zeros(3), zeros(3), "Ne"),
    Molecule(28.0134 * PHYSICS.u, 65e-12, zeros(3), zeros(3), "N2"),
    Molecule(31.9988 * PHYSICS.u, 60e-12, zeros(3), zeros(3), "O2"),
    Molecule(39.948 * PHYSICS.u, 1.88e-10, zeros(3), zeros(3), "Ar"),
]

function compute_next_pos!(mol, dt, g = PHYSICS.g)
    acc = [0, 0, -g]
    mol.velocity = mol.velocity + acc * dt
    mol.pos = mol.pos + mol.velocity * dt
end

function find_empty_pos(existing, mol::Molecule, domain::Domain)
    while true
        pos = ([domain.Lx, domain.Ly, domain.Lz] ./ 2 .- mol.radius) .* (2 .* rand(3) .- 1)
        if isempty(existing) || all(norm(m.pos .- pos) >= m.radius + mol.radius for m in existing)
            return pos
        end
    end
end

function generate_mol(n::Int, mol_types::Vector{Molecule}, domain::Domain; speed=1400.0)
    molecules = []
    for _ in 1:n
        new_mol = copy(rand(mol_types))
        v = 2 .* rand(3) .- 1
        new_mol.velocity = speed .* v ./ norm(v)
        new_mol.pos = find_empty_pos(molecules, new_mol, domain)
        push!(molecules, new_mol)
    end
    return molecules
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
    resolved = Set{Int}()
    for i in 1:length(molecules)
        for j in i+1:length(molecules)
            if !in(i, resolved) && !in(j, resolved) && are_colliding(molecules[i], molecules[j])
                resolve_collision!(molecules[i], molecules[j])
                push!(resolved, i)
                push!(resolved, j)
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
