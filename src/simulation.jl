const PLOT_SIZE = (800, 800)

const MOL_TYPES = [
    Molecule(4.0026 * PHYSICS.u, 1.1e-10, zeros(3), zeros(3), "He"),
    Molecule(20.1797 * PHYSICS.u, 38e-12, zeros(3), zeros(3), "Ne"),
    Molecule(28.0134 * PHYSICS.u, 65e-12, zeros(3), zeros(3), "N2"),
    Molecule(31.9988 * PHYSICS.u, 60e-12, zeros(3), zeros(3), "O2"),
]

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

function col_sim()
    dt = 1e-13
    duration = 40e-12
    domain = Domain(2e-9, 2e-9, 2e-9)
    lims = (-domain.Lx / 2, domain.Lx / 2)
    n_mols = 15

    t = 0.0
    molecules = generate_mol(n_mols, [MOL_TYPES[1]], domain)

    anim = @animate while t < duration
        p = plot3d(title="t = $(round(t*1e12, digits=2)) [ps]", xlims=lims, ylims=lims, zlims=lims, size=PLOT_SIZE)

        for mol in molecules
            scatter3d!(p, [mol.pos[1]], [mol.pos[2]], [mol.pos[3]], label=false)
            compute_next_pos!(mol, dt)
            reflect_walls!(mol, domain)
        end

        resolve_collisions!(molecules)

        t += dt
    end

    gif(anim, "export/col_sim.gif")
end

function single_col()
    dt = 1e-13
    duration = 10e-12
    domain = Domain(2e-9, 2e-9, 2e-9)
    lims = (-domain.Lx / 2, domain.Lx / 2)

    r_He = 31e-12
    m_He = MOL_TYPES[1].mass

    r_Ne = MOL_TYPES[2].radius
    m_Ne = MOL_TYPES[2].mass

    cases = [
        (
            Molecule(m_He, r_He, [-0.4e-9, 0.0, 0.0], [200.0, 0.0, 0.0], "He"),
            Molecule(m_He, r_He, [0.4e-9, 0.0, 0.0], [-200.0, 0.0, 0.0], "He"),
            "frontal shock, same masses",
            ),
        (
            Molecule(m_He, r_He, [-0.4e-9, 0.0, 0.0], [200.0, 0.0, 0.0], "He"),
            Molecule(m_Ne, r_Ne, [0.4e-9, 0.0, 0.0], [0.0, 0.0, 0.0], "Ne"),
            "frontal shock, different masses",
        ),
        (
            Molecule(m_He, r_He, [-0.4e-9, 0.0, 0.0], [200.0, 0.0, 0.0], "He"),
            Molecule(m_He, r_He, [0.0, (r_He + r_He) * sqrt(2) / 2, 0.0], [0.0, 0.0, 0.0], "He"),
            "oblique shock, same masses",
        ),
    ]

    for (i, (mol1, mol2, title)) in enumerate(cases)
        t = 0.0

        anim = @animate while t < duration
            p = plot(
                xlims=lims,
                ylims=lims,
                size=PLOT_SIZE,
            )

            scatter!(p,
                [mol1.pos[1]],
                [mol1.pos[2]],
                markersize=PLOT_SIZE[1] * 0.75 * mol1.radius / (domain.Lx / 2),
                color=:red,
                label=mol1.formula
            )

            scatter!(p,
                [mol2.pos[1]],
                [mol2.pos[2]],
                markersize=PLOT_SIZE[1] * 0.75 * mol2.radius / (domain.Lx / 2),
                color=:blue,
                label=mol2.formula
            )

            compute_next_pos!(mol1, dt)
            compute_next_pos!(mol2, dt)

            if are_colliding(mol1, mol2)
                resolve_collision!(mol1, mol2)
            end

            t += dt
        end

        gif(anim, "export/single_col$(i).gif")
    end
end

function he_sim()
    dt = 10e-14
    duration = 2e-11
    domain = Domain(1e-8, 1e-8, 1e-8)
    lims = (-domain.Lx / 2, domain.Lx / 2)
    n_mols = 400

    t = 0.0
    molecules = generate_mol(n_mols, [MOL_TYPES[1]], domain)

    v_means = []
    times = []
    velos_norms = []
    alphas = []
    betas = []

    anim = @animate while t < duration
        p = plot3d(title="t = $(round(t*1e12, digits=2)) [ps]", xlims=lims, ylims=lims, zlims=lims, size=PLOT_SIZE)

        velos = []
        velos_2 = []

        for mol in molecules
            scatter3d!(p, [mol.pos[1]], [mol.pos[2]], [mol.pos[3]], label=false)
            compute_next_pos!(mol, dt)
            reflect_walls!(mol, domain)
            push!(velos, norm(mol.velocity))
            push!(velos_2, norm(mol.velocity)^2)
        end

        velo_mean = mean(velos)
        velo2_mean = mean(velos_2)
        alpha = MOL_TYPES[1].mass * velo2_mean / (3 * PHYSICS.k_g)
        beta = n_mols * MOL_TYPES[1].mass * velo2_mean / (3 * volume(domain))

        push!(times, t)
        push!(v_means, velo_mean)
        push!(velos_norms, velos)
        push!(alphas, alpha)
        push!(betas, beta)

        resolve_collisions!(molecules)

        t += dt
    end

    # save the simulation animation
    gif(anim, "export/he_sim.gif")

    # display velocity mean during time simulation
    fig = plot(times, v_means, label=false, title="velocity mean", size=PLOT_SIZE)
    savefig(fig, "export/velo_mean.png")

    # display alpha during time simulation
    fig = plot(times, alphas, label=false, title="alpha (temperature)", ylims=(minimum(alphas) * 0.99, maximum(alphas) * 1.01), size=PLOT_SIZE)
    savefig(fig, "export/alpha.png")

    # display beta during time simulation
    fig = plot(times, betas, label=false, title="beta (pressure)", ylims=(minimum(betas) * 0.99, maximum(betas) * 1.01), size=PLOT_SIZE)
    savefig(fig, "export/beta.png")

    # display beginning distribution of mean during simulation
    fig = histogram(
        velos_norms[10],
        label=false,
        bins=40,
        title="10th velo distrib",
        xlims=(-500, 4000),
        ylims=(0, 200),
        size=PLOT_SIZE,
    )
    savefig(fig, "export/velo_10th.png")

    # display ending distribution of mean during simulation
    fig = histogram(
        last(velos_norms),
        label=false,
        bins=40,
        title="last velo distrib",
        xlims=(-500, 4000),
        ylims=(0, 200),
        size=PLOT_SIZE,
    )
    savefig(fig, "export/velo_last.png")

    anim = @animate for i in 1:length(velos_norms)
        histogram(
            velos_norms[i],
            label=false,
            bins=40,
            title="t [ps] = $(round(times[i]*1e12, digits=2))",
            xlims=(-500, 4000),
            ylim=(0, 120),
            size=PLOT_SIZE,
        )
    end
    gif(anim, "export/velo_distrib.gif")

    println()
end

function simulation()
    col_sim()
    single_col()
    he_sim()
end
