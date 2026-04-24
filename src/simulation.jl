const PLOT_SIZE = (800, 800)

function simple_collision()
    export_path = "export/simple_collision"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

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

    gif(anim, "$export_path/sim.gif")
end

function single_collision()
    export_path = "export/single_collision"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

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

        gif(anim, "$export_path/$(i).gif")
    end
end

function he()
    export_path = "export/he"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

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
    pos_z_all = []

    anim = @animate while t < duration
        p = plot3d(title="t = $(round(t*1e12, digits=2)) [ps]", xlims=lims, ylims=lims, zlims=lims, size=PLOT_SIZE)

        velos = []
        velos_2 = []
        pos_z = []

        for mol in molecules
            scatter3d!(p, [mol.pos[1]], [mol.pos[2]], [mol.pos[3]], label=false)
            compute_next_pos!(mol, dt)
            reflect_walls!(mol, domain)
            push!(velos, norm(mol.velocity))
            push!(velos_2, norm(mol.velocity)^2)
            push!(pos_z, mol.pos[3])
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
        push!(pos_z_all, pos_z)

        resolve_collisions!(molecules)

        t += dt
    end

    # save the simulation animation
    gif(anim, "$export_path/_sim.gif")

    # display velocity mean during time simulation
    fig = plot(times, v_means, label=false, title="velocity mean", size=PLOT_SIZE)
    savefig(fig, "$export_path/velo_mean.png")

    # display alpha during time simulation
    fig = plot(times, alphas, label=false, title="alpha (temperature)", ylims=(minimum(alphas) * 0.99, maximum(alphas) * 1.01), size=PLOT_SIZE)
    savefig(fig, "$export_path/alpha.png")

    # display beta during time simulation
    fig = plot(times, betas, label=false, title="beta (pressure)", ylims=(minimum(betas) * 0.99, maximum(betas) * 1.01), size=PLOT_SIZE)
    savefig(fig, "$export_path/beta.png")

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
    savefig(fig, "$export_path/velo_10th.png")

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
    savefig(fig, "$export_path/velo_last.png")

    anim = @animate for i in 1:length(velos_norms)
        histogram(
            velos_norms[i],
            label=false,
            bins=40,
            title="t [ps] = $(round(times[i]*1e12, digits=2))",
            xlims=(-500, 4000),
            ylims=(0, 120),
            size=PLOT_SIZE,
        )
    end
    gif(anim, "$export_path/velo_distrib.gif")

    anim = @animate for i in 1:length(pos_z_all)
        histogram(
            pos_z_all[i],
            label=false,
            bins=40,
            normalize=:probability,
            title="t [ps] = $(round(times[i]*1e12, digits=2))",
            xlims=lims,
            ylims=(0, 1)
        )
    end
    gif(anim, "$export_path/pos_z.gif")

    println()
end

function gravity_time(d, temp, mol)
    # define gravity velocity
    v = (mol.mass * PHYSICS.g * PHYSICS.D) / (PHYSICS.k_g * temp)

    # compute time in seconds
    s = d * 1 / v
    println("t [s]    : ", s)

    # convert to year
    h = s / 60 / 24 / 365
    println("t [year] : ", h)
end

function he_with_g()
    export_path = "export/he_with_g"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

    dt = 1e-13
    duration = 5e-11
    domain = Domain(2e-8, 2e-8, 2e-8)
    lims = (-domain.Lx / 2, domain.Lx / 2)
    n_mols = 1_000
    g = 9.81e13

    t = 0.0
    molecules = generate_mol(n_mols, [MOL_TYPES[1]], domain; speed=1367.0)

    pos_z_all = []
    times = []

    while t < duration
        pos_z = Float64[]

        for mol in molecules
            compute_next_pos!(mol, dt, g)
            reflect_walls!(mol, domain)
            push!(pos_z, mol.pos[3])
        end

        resolve_collisions!(molecules)

        push!(pos_z_all, pos_z)
        push!(times, t)

        t += dt
    end

    n_bins = 40

    # posz distribution from final state
    fig = histogram(
        last(pos_z_all),
        label=false,
        bins=n_bins,
        normalize=:probability,
        title="posz final distrib",
        xlabel="z [m]",
        ylabel="probability",
        size=PLOT_SIZE
    )
    savefig(fig, "$export_path/posz.png")

    # pressure profile from final state
    z_min = -domain.Lz / 2
    z_max = domain.Lz / 2
    dz = (z_max - z_min) / n_bins
    V_bin = domain.Lx * domain.Ly * dz
    mol_mass = MOL_TYPES[1].mass

    z_centers = [z_min + (i - 0.5) * dz for i in 1:n_bins]
    pressures = zeros(n_bins)

    for mol in molecules
        bin_idx = clamp(floor(Int, (mol.pos[3] - z_min) / dz) + 1, 1, n_bins)
        v2 = sum(mol.velocity .^ 2)
        pressures[bin_idx] += mol_mass * v2 / 3
    end
    pressures ./= V_bin

    fig = plot(
        z_centers,
        pressures,
        label=false,
        title="pressure z",
        xlabel="z [m]",
        ylabel="P [Pa]",
        size=PLOT_SIZE,
    )
    savefig(fig, "$export_path/pressure_z.png")

    println()
end

function he_ar()
    export_path = "export/he_ar"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

    dt = 1e-14
    duration = 5e-11
    domain = Domain(2e-8, 2e-8, 2e-8)
    lims = (-domain.Lx / 2, domain.Lx / 2)
    g = 9.81e13

    n_he = 400
    n_ar = 200

    he_type = MOL_TYPES[findfirst(m -> m.formula == "He", MOL_TYPES)]
    ar_type = MOL_TYPES[findfirst(m -> m.formula == "Ar", MOL_TYPES)]

    he_mols = generate_mol(n_he, [he_type], domain; speed=789.45)
    ar_mols = generate_mol(n_ar, [ar_type], domain; speed=249.88)
    molecules = [he_mols..., ar_mols...]

    t = 0.0
    step = 0
    anim_every = 30
    anim = Animation()

    while t < duration
        for mol in molecules
            compute_next_pos!(mol, dt, g)
            reflect_walls!(mol, domain)
        end
        resolve_collisions!(molecules)

        if step % anim_every == 0
            p = plot3d(
                title="t = $(round(t*1e12, digits=2)) [ps]",
                xlims=lims, ylims=lims, zlims=lims,
                size=PLOT_SIZE,
            )

            he_x = [m.pos[1] for m in molecules if m.formula == "He"]
            he_y = [m.pos[2] for m in molecules if m.formula == "He"]
            he_z = [m.pos[3] for m in molecules if m.formula == "He"]

            ar_x = [m.pos[1] for m in molecules if m.formula == "Ar"]
            ar_y = [m.pos[2] for m in molecules if m.formula == "Ar"]
            ar_z = [m.pos[3] for m in molecules if m.formula == "Ar"]

            scatter3d!(p, he_x, he_y, he_z, label="He", color=:blue, ms=2)
            scatter3d!(p, ar_x, ar_y, ar_z, label="Ar", color=:red, ms=3)

            frame(anim, p)
        end

        step += 1
        t += dt
    end

    # Q9.1
    gif(anim, "$export_path/_sim.gif")

    # Q9.2
    he_z_final = [m.pos[3] for m in molecules if m.formula == "He"]
    ar_z_final = [m.pos[3] for m in molecules if m.formula == "Ar"]

    fig = histogram(
        he_z_final,
        label="He",
        bins=40,
        normalize=:probability,
        color=:blue, alpha=0.5,
        title="z distribution per species (Maxwell–Boltzmann)",
        xlabel="z [m]",
        ylabel="probability",
        size=PLOT_SIZE,
    )
    histogram!(fig, ar_z_final,
        label="Ar",
        bins=40,
        normalize=:probability,
        color=:red, alpha=0.5,
    )
    savefig(fig, "$export_path/distrib.png")

    # Q9.4
    n_bins = 40
    z_min = -domain.Lz / 2
    z_max = domain.Lz / 2
    dz = (z_max - z_min) / n_bins
    V_bin = domain.Lx * domain.Ly * dz
    z_centers = [z_min + (i - 0.5) * dz for i in 1:n_bins]

    mv2_sum = zeros(n_bins)
    p_bin = zeros(n_bins)
    count_bin = zeros(Int, n_bins)

    for mol in molecules
        bin_idx = clamp(floor(Int, (mol.pos[3] - z_min) / dz) + 1, 1, n_bins)
        v2 = sum(mol.velocity .^ 2)
        mv2_sum[bin_idx] += mol.mass * v2
        p_bin[bin_idx] += mol.mass * v2 / 3
        count_bin[bin_idx] += 1
    end

    p_bin ./= V_bin
    mv2_mean = [count_bin[i] > 0 ? mv2_sum[i] / count_bin[i] : 0.0 for i in 1:n_bins]
    T_bin = mv2_mean ./ (3 * PHYSICS.k_g)

    fig = plot(z_centers, mv2_mean,
        label=false,
        title="mv2 and postion in z",
        xlabel="z [m]", ylabel="<mv²> [J]",
        ylims=(minimum(mv2_mean) * 0.99, maximum(mv2_mean) * 1.01),
        size=PLOT_SIZE,
    )
    savefig(fig, "$export_path/z-mv2.png")

    fig = plot(z_centers, p_bin,
        label=false,
        title="pressure and postion in z",
        xlabel="z [m]", ylabel="P [Pa]",
        size=PLOT_SIZE,
    )
    savefig(fig, "$export_path/z-pressure.png")

    fig = plot(z_centers, T_bin,
        label=false,
        title="temperature and postion in z",
        xlabel="z [m]", ylabel="T [K]",
        size=PLOT_SIZE,
    )
    savefig(fig, "$export_path/z-temp.png")

    he_v = [m.velocity for m in molecules if m.formula == "He"]
    ar_v = [m.velocity for m in molecules if m.formula == "Ar"]

    he_vnorm = [norm(v) for v in he_v]
    ar_vnorm = [norm(v) for v in ar_v]
    he_vx = [v[1] for v in he_v]
    he_vy = [v[2] for v in he_v]
    he_vz = [v[3] for v in he_v]
    ar_vx = [v[1] for v in ar_v]
    ar_vy = [v[2] for v in ar_v]
    ar_vz = [v[3] for v in ar_v]

    n_bins = 40
    function velo_subplot(he_data, ar_data, title_str, xlabel_str)
        p = histogram(he_data,
            label="He", bins=n_bins, normalize=:probability,
            color=:blue, alpha=0.5,
            title=title_str, xlabel=xlabel_str, ylabel="probability",
        )
        histogram!(p, ar_data,
            label="Ar", bins=n_bins, normalize=:probability,
            color=:red, alpha=0.5,
        )
        return p
    end

    p_norm = velo_subplot(he_vnorm, ar_vnorm, "||v|| distribution (Chi-squared)", "||v|| [m/s]")
    p_vx = velo_subplot(he_vx, ar_vx, "vx distribution", "vx [m/s]")
    p_vy = velo_subplot(he_vy, ar_vy, "vy distribution", "vy [m/s]")
    p_vz = velo_subplot(he_vz, ar_vz, "vz distribution", "vz [m/s]")

    fig = plot(p_norm, p_vx, p_vy, p_vz, layout=(2, 2), size=(PLOT_SIZE[1] * 2, PLOT_SIZE[2] * 2))
    savefig(fig, "$export_path/velo_distrib.png")

    println()
end

function stability()
    export_path = "export/stability"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

    dt = 10e-14
    duration = 4e-10
    domain = Domain(1e-8, 1e-8, 1e-8)
    n_mols = 400

    t = 0.0
    molecules = generate_mol(n_mols, [MOL_TYPES[1]], domain)

    speed = 1400.0
    for mol in molecules
        mol.velocity = [rand([-1.0, 1.0]) * speed, 0.0, 0.0]
    end

    times = Float64[]
    vx2_mean = Float64[]
    vy2_mean = Float64[]
    vz2_mean = Float64[]

    while t < duration
        vx2 = Float64[]
        vy2 = Float64[]
        vz2 = Float64[]
        for mol in molecules
            compute_next_pos!(mol, dt)
            reflect_walls!(mol, domain)
            push!(vx2, mol.velocity[1]^2)
            push!(vy2, mol.velocity[2]^2)
            push!(vz2, mol.velocity[3]^2)
        end
        resolve_collisions!(molecules)

        push!(times, t)
        push!(vx2_mean, mean(vx2))
        push!(vy2_mean, mean(vy2))
        push!(vz2_mean, mean(vz2))

        t += dt
    end

    window = 100
    treshold = 0.05

    function rel_drift(values, W)
        n = length(values)
        out = fill(NaN, n)
        for i in 2W:n
            μ_prev = mean(values[i-2W+1:i-W])
            μ_curr = mean(values[i-W+1:i])
            out[i] = abs(μ_curr - μ_prev) / abs(μ_curr)
        end
        return out
    end

    drift_x = rel_drift(vx2_mean, window)
    drift_y = rel_drift(vy2_mean, window)
    drift_z = rel_drift(vz2_mean, window)

    fig = plot(times, vx2_mean, label="<vx²>", title="velocities repartition",
        xlabel="t [s]", ylabel="<v²> [m²/s²]", size=PLOT_SIZE)
    plot!(fig, times, vy2_mean, label="<vy²>")
    plot!(fig, times, vz2_mean, label="<vz²>")
    savefig(fig, "$export_path/repartition.png")

    fig = plot(times, drift_x, label="x", title="stability factor : sliding window",
        xlabel="t [s]", ylabel="|Δμ|/μ", size=PLOT_SIZE)
    plot!(fig, times, drift_y, label="y")
    plot!(fig, times, drift_z, label="z")
    hline!(fig, [treshold], label="seuil $(treshold)", linestyle=:dash)
    savefig(fig, "$export_path/stability.png")

    println()
end

function launch_simulation()
    # simple_collision()
    # single_collision()

    # he()
    # gravity_time(1, 26.85, MOL_TYPES[1])
    # he_with_g()

    he_ar()

    stability()
end
