const DT = 1e-13
const DOMAIN = Domain(2e-9, 2e-9, 2e-9)
const N_MOLS = 15
const DURATION = 40e-12

const LIMS = (-DOMAIN.Lx / 2, DOMAIN.Lx / 2)
const PLOT_SIZE = (800, 800)

const MOL_TYPES = [
    Molecule(4.0026 * PHYSICS.u, 1.1e-10, zeros(3), zeros(3), "He"),
    Molecule(20.1797 * PHYSICS.u, 38e-12, zeros(3), zeros(3), "Ne"),
    Molecule(28.0134 * PHYSICS.u, 65e-12, zeros(3), zeros(3), "N2"),
    Molecule(31.9988 * PHYSICS.u, 60e-12, zeros(3), zeros(3), "O2"),
]

function find_empty_pos(existing, mol::Molecule)
    while true
        pos = ([DOMAIN.Lx, DOMAIN.Ly, DOMAIN.Lz] ./ 2 .- mol.radius) .* (2 .* rand(3) .- 1)
        if isempty(existing) || all(norm(m.pos .- pos) >= m.radius + mol.radius for m in existing)
            return pos
        end
    end
end

function generate_mol(n::Int)
    molecules = []
    for _ in 1:n
        new_mol = copy(MOL_TYPES[1])
        # new_mol.velocity = 100.0 .* (2 .* rand(3) .- 1)
        v = 2 .* rand(3) .- 1
        new_mol.velocity = 1400.0 .* v ./ norm(v)
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
            reflect_walls!(mol, DOMAIN)
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
                markersize=PLOT_SIZE[1] * 0.75 * mol1.radius / (DOMAIN.Lx / 2),
                color=:red,
                label=mol1.formula
            )

            scatter!(p,
                [mol2.pos[1]],
                [mol2.pos[2]],
                markersize=PLOT_SIZE[1] * 0.75 * mol2.radius / (DOMAIN.Lx / 2),
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

function single_gas()
    domain = Domain(5e-9, 5e-9, 5e-9)
    lims = (-domain.Lx / 2, domain.Lx / 2)
    n_mols = 400
    v0 = 1400
    dt = 10e-14
    duration = 2e-11

    t = 0.0
    molecules = generate_mol(n_mols)

    v_means = []
    times = []
    velos_norms = []
    alphas = []
    betas = []

    anim = @animate while t < duration
        p = plot3d(
            title="t = $(round(t*1e12, digits=2)) [ps]",
            xlims=lims,
            ylims=lims,
            zlims=lims,
            size=PLOT_SIZE
        )

        velos = []
        for mol in molecules
            scatter3d!(p, [mol.pos[1]], [mol.pos[2]], [mol.pos[3]], label=false)
            compute_next_pos!(mol, dt)
            reflect_walls!(mol, domain)
            push!(velos, norm(mol.velocity))
        end
        velo_mean = mean(velos)
        push!(v_means, velo_mean)
        push!(velos_norms, velos)
        push!(times, t)

        alpha = (MOL_TYPES[1].mass * velo_mean^2) / (3 * PHYSICS.k_g)
        push!(alphas, alpha)

        beta = (n_mols * MOL_TYPES[1].mass * velo_mean^2) / (3 * volume(domain))
        push!(betas, beta)

        for i in 1:length(molecules)
            for j in i+1:length(molecules)
                if are_colliding(molecules[i], molecules[j])
                    resolve_collision!(molecules[i], molecules[j])
                end
            end
        end

        t += DT
    end

    # save the simulation animation
    gif(anim, "export/single_gas_sim.gif")

    # display velocity mean during time simulation
    fig = plot(times, v_means, label=false, title="velocity mean", size=(800, 800))
    display(fig)

    # display alpha during time simulation
    fig = plot(times, alphas, label=false, title="alphas", size=(800, 800))
    display(fig)

    # display beta during time simulation
    fig = plot(times, betas, label=false, title="betas", size=(800, 800))
    display(fig)

    # display beginning distribution of mean during simulation
    fig = histogram(
        velos_norms[3],
        label=false,
        bins=40,
        title="beginning velo distribution",
        xlims=(-500, 4000),
        ylims=(0, 200),
        size=(800, 800),
    )
    display(fig)

    # display ending distribution of mean during simulation
    fig = histogram(
        last(velos_norms),
        label=false,
        bins=40,
        title="ending velo distribution",
        xlims=(-500, 4000),
        ylims=(0, 200),
        size=(800, 800),
    )
    display(fig)

    anim = @animate for i in 1:length(velos_norms)
        histogram(
            velos_norms[i],
            label=false,
            bins=40,
            title="t [ps] = $(round(times[i]*1e12, digits=2))",
            xlims=(-500, 4000),
            ylim=(0, 120),
            size=(800, 800),
        )
    end

    gif(anim, "export/velocity_distribution.gif")
end

function simulation()
    # simple_sim()
    # collision_sim()
    # single_collision()
    single_gas()
end
