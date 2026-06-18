# Q13.2 — Réémission thermique selon la distribution de Maxwell-Boltzmann
# vx, vy ~ N(0, σ) — composantes tangentielles (distribution d'équilibre)
# vz ~ Rayleigh   — composante normale (flux sortant ∝ vz * f(vz))
# σ = √(kB * T / m)
function thermal_velocity(mass, T)
    σ = sqrt(PHYSICS.k_g * T / mass)
    vx = randn() * σ
    vy = randn() * σ
    vz = σ * sqrt(-2 * log(rand()))   # distribution de Rayleigh
    return (vx, vy, vz)
end

# Q13.2 — Condition de bord mixte :
#   x, y : réflexion spéculaire (comme reflect_walls!)
#   z+   : paroi chaude à T_top — réémission thermique vers l'intérieur (vz < 0)
#   z-   : paroi froide à T_bot — réémission thermique vers l'intérieur (vz > 0)
function reflect_thermal_walls!(mol::Molecule, domain::Domain, T_top, T_bot)
    half = [domain.Lx, domain.Ly, domain.Lz] ./ 2

    # réflexion spéculaire pour x et y
    for i in 1:2
        lo = -half[i] + mol.radius
        hi =  half[i] - mol.radius
        if mol.pos[i] < lo
            mol.pos[i] = 2 * lo - mol.pos[i]
            mol.velocity[i] = abs(mol.velocity[i])
        elseif mol.pos[i] > hi
            mol.pos[i] = 2 * hi - mol.pos[i]
            mol.velocity[i] = -abs(mol.velocity[i])
        end
    end

    # paroi z- (bord bas, T_bot) → réémission vers +z
    lo_z = -half[3] + mol.radius
    if mol.pos[3] < lo_z
        mol.pos[3] = 2 * lo_z - mol.pos[3]
        vx, vy, vz = thermal_velocity(mol.mass, T_bot)
        mol.velocity = [vx, vy, abs(vz)]   # vz > 0 : vers l'intérieur
    end

    # paroi z+ (bord haut, T_top) → réémission vers -z
    hi_z = half[3] - mol.radius
    if mol.pos[3] > hi_z
        mol.pos[3] = 2 * hi_z - mol.pos[3]
        vx, vy, vz = thermal_velocity(mol.mass, T_top)
        mol.velocity = [vx, vy, -abs(vz)]  # vz < 0 : vers l'intérieur
    end
end
