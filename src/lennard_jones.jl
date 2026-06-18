# Q14.3 — Force de Lennard-Jones exercée sur mol_i par tous les j ≠ i (coupure à 3σ)
# F_i = Σ_j 24ε/r² * [(σ/r)^6 - 2(σ/r)^12] * r_ij
# avec r_ij = r_j - r_i  (de i vers j)
# → répulsif pour r < 2^(1/6)σ, attractif pour r > 2^(1/6)σ
function lj_total_force(mol_i::Molecule, molecules, sigma::Float64, epsilon::Float64)
    F = zeros(3)
    cutoff = 3 * sigma
    for mol_j in molecules
        mol_j === mol_i && continue
        r_vec = mol_j.pos .- mol_i.pos
        r = norm(r_vec[1:2])   # distance 2D (simulation plane xy)
        (r < 1e-15 || r > cutoff) && continue
        sr6 = (sigma / r)^6
        coeff = 24 * epsilon / r^2 * (sr6 - 2 * sr6^2)
        F[1] += coeff * r_vec[1]
        F[2] += coeff * r_vec[2]
    end
    return F
end

# Q14.4 — Température d'un système 2D : T = m<v²> / (2 kB)
# (facteur 2 car 2 degrés de liberté translationels, vs 3 en 3D)
# Q14.5 — Retourne T_cal
function temperature_2d(molecules)
    v2_mean = mean(m.velocity[1]^2 + m.velocity[2]^2 for m in molecules)
    return molecules[1].mass * v2_mean / (2 * PHYSICS.k_g)
end

# Q14.6 — Rescaling des vitesses : vi ← vi * √(T_ref / T_cal)
function rescale_velocities!(molecules, T_ref)
    T_cal = temperature_2d(molecules)
    T_cal < 1e-30 && return
    factor = sqrt(T_ref / T_cal)
    for mol in molecules
        mol.velocity[1] *= factor
        mol.velocity[2] *= factor
    end
end

# Initialisation Maxwell-Boltzmann 2D : vx,vy ~ N(0, √(kB*T/m)), vz = 0
function maxwell_boltzmann_2d(mass::Float64, T::Float64)
    σ = sqrt(PHYSICS.k_g * T / mass)
    return [randn() * σ, randn() * σ, 0.0]
end

# Position aléatoire dans le plan xy sans chevauchement (distance ≥ 0.9σ_lj entre particules)
function find_empty_pos_2d(existing, Lx, Ly, sigma_lj)
    while true
        x = Lx / 2 * (2 * rand() - 1)
        y = Ly / 2 * (2 * rand() - 1)
        pos = [x, y, 0.0]
        ok = all(norm(m.pos[1:2] .- pos[1:2]) >= 0.9 * sigma_lj for m in existing)
        (isempty(existing) || ok) && return pos
    end
end

# Génère n molécules 2D LJ avec initialisation MB et sans chevauchement
function generate_mol_2d_lj(n::Int, mol_type::Molecule, Lx::Float64, Ly::Float64, sigma_lj::Float64, T_ref::Float64)
    molecules = Molecule[]
    for _ in 1:n
        mol = copy(mol_type)
        mol.velocity = maxwell_boltzmann_2d(mol.mass, T_ref)
        mol.pos = find_empty_pos_2d(molecules, Lx, Ly, sigma_lj)
        push!(molecules, mol)
    end
    return molecules
end

# Réflexion spéculaire 2D (pour LJ, pas de rayon dur — réflexion au bord du domaine)
function reflect_walls_2d!(mol::Molecule, Lx::Float64, Ly::Float64)
    half_x = Lx / 2
    half_y = Ly / 2
    if mol.pos[1] < -half_x
        mol.pos[1] = -2 * half_x - mol.pos[1]
        mol.velocity[1] = abs(mol.velocity[1])
    elseif mol.pos[1] > half_x
        mol.pos[1] = 2 * half_x - mol.pos[1]
        mol.velocity[1] = -abs(mol.velocity[1])
    end
    if mol.pos[2] < -half_y
        mol.pos[2] = -2 * half_y - mol.pos[2]
        mol.velocity[2] = abs(mol.velocity[2])
    elseif mol.pos[2] > half_y
        mol.pos[2] = 2 * half_y - mol.pos[2]
        mol.velocity[2] = -abs(mol.velocity[2])
    end
end

# Q15.1 — Profil de température cible pour la condensation :
#   - N1 premiers pas : T = T1 (phase gazeuse)
#   - pas N1+1 à N1+N2 : T décroît linéairement de T1 à T2
#   - après N1+N2 : T = T2 (phase solide)
function target_temperature(step::Int, T1::Float64, N1::Int, T2::Float64, N2::Int)
    if step <= N1
        return T1
    elseif step <= N1 + N2
        frac = (step - N1) / N2
        return T1 + frac * (T2 - T1)
    else
        return T2
    end
end
