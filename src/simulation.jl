const DT = 1e-13
const BOX = 1e-9
const N_MOLS = 15
const DURATION = 40e-12

const LIMS = (-BOX, BOX)
const PLOT_SIZE = (800, 800)

const MOL_TYPES = [
    Molecule(4.0026 * PHYSICS.u, 31e-12, zeros(3), zeros(3), "He"),
    Molecule(20.1797 * PHYSICS.u, 38e-12, zeros(3), zeros(3), "Ne"),
    Molecule(28.0134 * PHYSICS.u, 65e-12, zeros(3), zeros(3), "N2"),
    Molecule(31.9988 * PHYSICS.u, 60e-12, zeros(3), zeros(3), "O2"),
]

function find_empty_pos(existing, mol::Molecule)
    while true
        pos = (BOX - mol.radius) .* (2 .* rand(3) .- 1)
        if isempty(existing) || all(norm(m.pos .- pos) >= m.radius + mol.radius for m in existing)
            return pos
        end
    end
end

function generate_mol(n::Int)
    molecules = []
    for _ in 1:n
        # new_mol = copy(MOL_TYPES[rand(1:length(MOL_TYPES))])
        new_mol = copy(MOL_TYPES[1])
        new_mol.velocity = 100.0 .* (2 .* rand(3) .- 1)
        new_mol.pos = find_empty_pos(molecules, new_mol)
        push!(molecules, new_mol)
    end
    return molecules
end

function simple_sim()
    t = 0.0
    molecules = generate_mol(N_MOLS)

    anim = @animate while t < DURATION
        p = plot3d(
            title="t = $(round(t*1e12, digits=2)) [ps]",
            xlims=LIMS,
            ylims=LIMS,
            zlims=LIMS,
            size=PLOT_SIZE
        )
        for mol in molecules
            scatter3d!(p, [mol.pos[1]], [mol.pos[2]], [mol.pos[3]], label=false)
            compute_next_pos!(mol, DT)
            reflect_walls!(mol, BOX)
        end
        t += DT
    end

    gif(anim, "export/simple_sim.gif")
end

function collision_sim()
    t = 0.0
    molecules = generate_mol(N_MOLS)

    anim = @animate while t < DURATION
        p = plot3d(
            title="t = $(round(t*1e12, digits=2)) [ps]",
            xlims=LIMS,
            ylims=LIMS,
            zlims=LIMS,
            size=PLOT_SIZE
        )

        for mol in molecules
            scatter3d!(p, [mol.pos[1]], [mol.pos[2]], [mol.pos[3]], label=false)
            compute_next_pos!(mol, DT)
            reflect_walls!(mol, BOX)
        end

        for i in 1:length(molecules)
            for j in i+1:length(molecules)
                if are_colliding(molecules[i], molecules[j])
                    resolve_collision!(molecules[i], molecules[j])
                end
            end
        end

        t += DT
    end

    gif(anim, "export/collision_sim.gif")
end

function single_collision()
    r_He, r_Ne = 31e-12, 38e-12
    m_He = 4.0026 * PHYSICS.u
    m_Ne = 20.1797 * PHYSICS.u

    cases = [
        (
            Molecule(m_He, r_He, [-0.4e-9, 0.0, 0.0], [100.0, 0.0, 0.0], "He"),
            Molecule(m_He, r_He, [0.4e-9, 0.0, 0.0], [-100.0, 0.0, 0.0], "He"),
            "cas 1 : frontal, masses égales (He+He)",
        ),
        (
            Molecule(m_He, r_He, [-0.4e-9, 0.0, 0.0], [100.0, 0.0, 0.0], "He"),
            Molecule(m_Ne, r_Ne, [0.4e-9, 0.0, 0.0], [0.0, 0.0, 0.0], "Ne"),
            "cas 2 : frontal, masses différentes (He+Ne)",
        ),
        (
            Molecule(m_He, r_He, [-0.4e-9, 0.0, 0.0], [100.0, 0.0, 0.0], "He"),
            Molecule(m_He, r_He, [0.0, (r_He + r_He) * sqrt(2) / 2, 0.0], [0.0, 0.0, 0.0], "He"),
            "cas 3 : oblique 45° (He+He)",
        ),
    ]

    duration = 10e-12

    for (i, (mol1, mol2, title)) in enumerate(cases)
        t = 0.0

        anim = @animate while t < duration
            p = plot(
                xlims=LIMS,
                ylims=LIMS,
                size=PLOT_SIZE,
            )

            scatter!(p,
                [mol1.pos[1]],
                [mol1.pos[2]],
                markersize=PLOT_SIZE[1] * 0.75 * mol1.radius / BOX,
                color=:red,
                label=mol1.formula
            )

            scatter!(p,
                [mol2.pos[1]],
                [mol2.pos[2]],
                markersize=PLOT_SIZE[1] * 0.75 * mol2.radius / BOX,
                color=:blue,
                label=mol2.formula
            )

            compute_next_pos!(mol1, DT)
            compute_next_pos!(mol2, DT)

            if are_colliding(mol1, mol2)
                resolve_collision!(mol1, mol2)
            end

            t += DT
        end

        gif(anim, "export/single_collision$(i).gif")
    end
end

function simulation()
    simple_sim()
    collision_sim()
    single_collision()
end
