const PLOT_SIZE = (800, 800)


# Q4.2 : cas de test unitaires pour la collision élastique
function single_collision()
    export_path = "export/q04"

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

        gif(anim, "$export_path/single_$(i).gif")
    end
end


# Q4.4 : simulation multi-molécules pour valider l'algorithme de collision
function simple_collision()
    export_path = "export/q04"
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

    gif(anim, "$export_path/multi_sim.gif")
end


# Q7.1 – Q7.5 : première simulation He — animation, vitesse moyenne, α (T), β (P)
function he()
    export_path = "export/q07"
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

    # Q7.1 — animation de la simulation
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

    gif(anim, "$export_path/_sim.gif")

    # Q7.2 — évolution de la vitesse moyenne
    fig = plot(times, v_means, label=false, title="velocity mean", size=PLOT_SIZE)
    savefig(fig, "$export_path/velo_mean.png")

    # Q7.4 — α = m<v²>/(3kB) → température
    fig = plot(times, alphas, label=false, title="alpha (temperature)", ylims=(minimum(alphas) * 0.99, maximum(alphas) * 1.01), size=PLOT_SIZE)
    savefig(fig, "$export_path/alpha.png")

    # Q7.5 — β = Nm<v²>/(3V) → pression
    fig = plot(times, betas, label=false, title="beta (pressure)", ylims=(minimum(betas) * 0.99, maximum(betas) * 1.01), size=PLOT_SIZE)
    savefig(fig, "$export_path/beta.png")

    # Q7.3 — distribution de la magnitude des vitesses au 10e pas
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

    # Q7.3 — distribution de la magnitude des vitesses au dernier pas
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
end


# Q8.3 : calcul du temps de parcours d'une molécule sous gravité seule
function gravity_time(d=1, temp=26.85, mol=MOL_TYPES[1])
    export_path = "export/q08"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

    # define gravity velocity
    v = (mol.mass * PHYSICS.g * PHYSICS.D) / (PHYSICS.k_g * temp)

    # compute time in seconds
    s = d * 1 / v

    # convert to year
    h = s / 60 / 24 / 365

    # save res to txt
    open("$export_path/gravity_time.txt", "w") do io
        println(io, "t [s]    : ", s)
        println(io, "t [year] : ", h)
    end
end


# Q8.4 : simulation He avec gravité exagérée — distribution z et pression
# Q8.6 : profil de pression en fonction de z
function he_with_g()
    export_path = "export/q08"

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

    # Q8.4 — distribution de probabilité de présence en z (état final)
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

    # Q8.6 — profil de pression en fonction de z
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
end


# Q9.1 : simulation He+Ar avec gravité — animation 3D différenciant les espèces
# Q9.2 : distribution de présence en z par espèce
# Q9.4 : <mv²>, pression et température en fonction de z
function he_ar()
    export_path = "export/q09"
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

    # Q9.1 — animation différenciant He (bleu) et Ar (rouge)
    gif(anim, "$export_path/_sim.gif")

    # Q9.2 — distribution de présence en z par espèce
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

    # Q9.4 — <mv²>, pression et température en fonction de z
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

    p_norm = velo_subplot(he_vnorm, ar_vnorm, "||v|| distribution (Chi)", "||v|| [m/s]")
    p_vx = velo_subplot(he_vx, ar_vx, "vx distribution", "vx [m/s]")
    p_vy = velo_subplot(he_vy, ar_vy, "vy distribution", "vy [m/s]")
    p_vz = velo_subplot(he_vz, ar_vz, "vz distribution", "vz [m/s]")

    fig = plot(p_norm, p_vx, p_vy, p_vz, layout=(2, 2), size=(PLOT_SIZE[1] * 2, PLOT_SIZE[2] * 2))
    savefig(fig, "$export_path/velo_distrib.png")
end


# Q10.2 : critère de stabilité — sliding window sur <vx²>, <vy²>, <vz²>
function stability()
    export_path = "export/q10"
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
            μ_prev = mean(values[(i-2W+1):(i-W)])
            μ_curr = mean(values[(i-W+1):i])
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
end


# Q12.2 : He 400 atomes concentrés dans 1/8 du domaine — entropie de Shannon H(t)
function shannon_half_domain()
    export_path = "export/q12"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

    dt = 1e-14
    duration = 5e-11
    Lxyz = 1e-8
    domain = Domain(Lxyz, Lxyz, Lxyz)
    n_mols = 400
    speed = 1400.0
    mol_type = MOL_TYPES[1]   # He

    # positions initiales dans [-5e-9, 0]³ (1/8 du domaine)
    molecules = Molecule[]
    while length(molecules) < n_mols
        mol = copy(mol_type)
        v = 2 .* rand(3) .- 1
        mol.velocity = speed .* v ./ norm(v)
        mol.pos = [-Lxyz / 2 + rand() * Lxyz / 2,
            -Lxyz / 2 + rand() * Lxyz / 2,
            -Lxyz / 2 + rand() * Lxyz / 2]
        push!(molecules, mol)
    end

    t = 0.0
    times = Float64[]
    entropies = Float64[]

    while t < duration
        for mol in molecules
            compute_next_pos!(mol, dt, 0.0)   # sans gravité
            reflect_walls!(mol, domain)
        end
        resolve_collisions!(molecules)

        push!(times, t)
        push!(entropies, total_entropy(molecules, domain))

        t += dt
    end

    fig = plot(times, entropies,
        label=false, title="Entropie de Shannon H(X,Y,Z,<v²>)",
        xlabel="t [s]", ylabel="H [bits]", size=PLOT_SIZE)
    savefig(fig, "$export_path/entropy_half.png")
    return molecules, domain, times, entropies
end


# Q12.3 : ouverture d'une paroi — domaine étendu sur x, nouvelle croissance de H
function shannon_open_wall()
    export_path = "export/q12"
    isdir(export_path) || mkpath(export_path)

    molecules, domain_small, times1, entropies1 = shannon_half_domain()

    dt = 1e-14
    n_extra = 5000
    domain_big = Domain(2 * domain_small.Lx, domain_small.Ly, domain_small.Lz)

    times2 = Float64[]
    entropies2 = Float64[]
    t = times1[end] + dt

    for _ in 1:n_extra
        for mol in molecules
            compute_next_pos!(mol, dt, 0.0)
            reflect_walls!(mol, domain_big)
        end
        resolve_collisions!(molecules)

        push!(times2, t)
        push!(entropies2, total_entropy(molecules, domain_big))

        t += dt
    end

    all_times = vcat(times1, times2)
    all_entropies = vcat(entropies1, entropies2)
    t_open = times1[end]

    fig = plot(all_times, all_entropies,
        label=false, title="Entropie — ouverture de paroi",
        xlabel="t [s]", ylabel="H [bits]", size=PLOT_SIZE)
    vline!(fig, [t_open], label="ouverture paroi", linestyle=:dash, color=:red)
    savefig(fig, "$export_path/entropy_open_wall.png")
end


# Q13.3 : He 500 atomes, parois z thermiques (T_top=700K, T_bot=300K)
# Q13.4 : évolution de l'entropie au cours du temps
# Q13.5 : profil de température T(z) à la fin de la simulation
function temp_gradient_sim()
    export_path = "export/q13"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

    dt = 1e-14
    duration = 1e-10
    domain = Domain(4e-9, 4e-9, 1e-8)
    n_mols = 500
    T_top = 700.0
    T_bot = 300.0

    t = 0.0
    molecules = generate_mol(n_mols, [MOL_TYPES[1]], domain; speed=1400.0)

    times = Float64[]
    entropies = Float64[]

    while t < duration
        for mol in molecules
            compute_next_pos!(mol, dt, 0.0)            # sans gravité
            reflect_thermal_walls!(mol, domain, T_top, T_bot)
        end
        resolve_collisions!(molecules)

        push!(times, t)
        push!(entropies, total_entropy(molecules, domain))

        t += dt
    end

    # Q13.4 — entropie au cours du temps
    fig = plot(times, entropies,
        label=false, title="Entropie — gradient de température",
        xlabel="t [s]", ylabel="H [bits]", size=PLOT_SIZE)
    savefig(fig, "$export_path/entropy.png")

    # Q13.5 — profil T(z) final
    n_bins = 20
    z_min = -domain.Lz / 2
    z_max = domain.Lz / 2
    dz = (z_max - z_min) / n_bins

    count_bin = zeros(Int, n_bins)
    mv2_bin = zeros(Float64, n_bins)

    for mol in molecules
        idx = clamp(floor(Int, (mol.pos[3] - z_min) / dz) + 1, 1, n_bins)
        v2 = mol.velocity[1]^2 + mol.velocity[2]^2 + mol.velocity[3]^2
        mv2_bin[idx] += mol.mass * v2
        count_bin[idx] += 1
    end

    z_centers = [z_min + (i - 0.5) * dz for i in 1:n_bins]
    T_profile = [count_bin[i] > 0 ? mv2_bin[i] / count_bin[i] / (3 * PHYSICS.k_g) : 0.0
                 for i in 1:n_bins]

    fig = plot(z_centers, T_profile,
        label=false, title="Profil de température T(z)",
        xlabel="z [m]", ylabel="T [K]", size=PLOT_SIZE)
    hline!(fig, [T_bot], label="T_bot=$(T_bot)K", linestyle=:dash, color=:blue)
    hline!(fig, [T_top], label="T_top=$(T_top)K", linestyle=:dash, color=:red)
    savefig(fig, "$export_path/T_vs_z.png")
end


# Q14.7 : boucle principale LJ — forces → Euler → réflexion → velocity rescaling
# Q14.8 : animation 2D (markersize=13)
# Q14.9 : comparer T_ref=10K (solide) et T_ref=40K (gaz/liquide)
function neon_lj_sim(T_ref=10.0, export_path="export/q14")
    label = string(Int(round(T_ref)))

    sigma = 2.74e-10
    epsilon = 4.91511044e-22
    dt = 1e-15
    duration = 5e-11
    Lx = Ly = 1e-8
    n_mols = 100

    ne_type = MOL_TYPES[findfirst(m -> m.formula == "Ne", MOL_TYPES)]
    molecules = generate_mol_2d_lj(n_mols, ne_type, Lx, Ly, sigma, T_ref)

    t = 0.0
    step = 0
    anim_every = max(1, floor(Int, duration / dt / 200))
    anim = Animation()
    lims = (-Lx / 2, Lx / 2)

    while t < duration
        # Q14.7.1 — calcul des forces LJ pour toutes les molécules
        forces = [lj_total_force(mol, molecules, sigma, epsilon) for mol in molecules]

        # Q14.7.2 — intégration Euler explicite
        for (mol, F) in zip(molecules, forces)
            mol.velocity[1] += F[1] / mol.mass * dt
            mol.velocity[2] += F[2] / mol.mass * dt
            mol.pos[1] += mol.velocity[1] * dt
            mol.pos[2] += mol.velocity[2] * dt
        end

        # Q14.7.3 — condition de bord (réflexion spéculaire 2D)
        for mol in molecules
            reflect_walls_2d!(mol, Lx, Ly)
        end

        # Q14.7.4 — contrôle de température (velocity rescaling)
        rescale_velocities!(molecules, T_ref)

        # Q14.8 — frame d'animation
        if step % anim_every == 0
            p = scatter(
                [m.pos[1] for m in molecules],
                [m.pos[2] for m in molecules],
                markersize=13, label=false,
                title="Ne LJ — T_ref=$(T_ref) K — t=$(round(t*1e12, digits=1)) ps",
                xlims=lims, ylims=lims, size=PLOT_SIZE,
            )
            frame(anim, p)
        end

        step += 1
        t += dt
    end

    gif(anim, "$export_path/T$(label).gif")
end


# Q15.2 : simulation de condensation solide Ne LJ (T décroît de 40K à 10K)
# Q15.3 : animation sur les 3 phases (gazeuse, transition, solide)
# Q15.4 : distribution spatiale de probabilité pour chaque phase
function condensation_sim()
    export_path = "export/q15"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)

    sigma = 2.74e-10
    epsilon = 4.91511044e-22
    dt = 1e-15
    Lx = Ly = 5e-9    # boîte resserrée (vs 1e-8) pour densifier et favoriser la nucléation
    n_mols = 100
    T1 = 40.0
    T2 = 10.0
    N1 = 15_000
    N2 = 15_000
    N_total = N1 + N2 + N2    # N2 pas supplémentaires à T2 après la transition

    ne_type = MOL_TYPES[findfirst(m -> m.formula == "Ne", MOL_TYPES)]
    molecules = generate_mol_2d_lj(n_mols, ne_type, Lx, Ly, sigma, T1)

    anim_every = max(1, N_total ÷ 200)
    anim = Animation()
    lims = (-Lx / 2, Lx / 2)

    # snapshots pour Q15.4 (fin de chaque phase)
    snap_gas = nothing
    snap_trans = nothing
    snap_solid = nothing

    for step in 1:N_total
        # Q15.1 — température cible selon la phase
        T_ref = target_temperature(step, T1, N1, T2, N2)

        forces = [lj_total_force(mol, molecules, sigma, epsilon) for mol in molecules]

        for (mol, F) in zip(molecules, forces)
            mol.velocity[1] += F[1] / mol.mass * dt
            mol.velocity[2] += F[2] / mol.mass * dt
            mol.pos[1] += mol.velocity[1] * dt
            mol.pos[2] += mol.velocity[2] * dt
        end

        for mol in molecules
            reflect_walls_2d!(mol, Lx, Ly)
        end

        rescale_velocities!(molecules, T_ref)

        # Q15.3 — frames d'animation
        if step % anim_every == 0
            T_cur = round(T_ref, digits=1)
            phase = step <= N1 ? "gaz" : (step <= N1 + N2 ? "transition" : "solide")
            p = scatter(
                [m.pos[1] for m in molecules],
                [m.pos[2] for m in molecules],
                markersize=13, label=false,
                title="Condensation Ne — $phase — T=$(T_cur) K — step=$step",
                xlims=lims, ylims=lims, size=PLOT_SIZE,
            )
            frame(anim, p)
        end

        # Q15.4 — capture des positions à la fin de chaque phase
        if step == N1
            snap_gas = [(m.pos[1], m.pos[2]) for m in molecules]
        elseif step == N1 + N2
            snap_trans = [(m.pos[1], m.pos[2]) for m in molecules]
        elseif step == N_total
            snap_solid = [(m.pos[1], m.pos[2]) for m in molecules]
        end
    end

    # Q15.3 — animation complète
    gif(anim, "$export_path/sim.gif")

    # Q15.4 — distributions spatiales par phase
    edges = range(-Lx / 2, Lx / 2, length=21)   # bins calés sur le domaine (échelle commune aux 3 phases)
    function plot_spatial_dist(snap, title_str, filename)
        xs = first.(snap)
        ys = last.(snap)
        p = histogram2d(xs, ys,
            bins=(edges, edges),
            color=:viridis,                  # colormap perceptuelle : contraste net amas/vide
            normalize=:probability,
            clims=(0.0, 0.06),               # plafond fixe → un amas dense ressort en jaune
            title=title_str, xlabel="x [m]", ylabel="y [m]",
            xlims=lims, ylims=lims, size=PLOT_SIZE)
        savefig(p, "$export_path/$filename")
    end

    snap_gas !== nothing && plot_spatial_dist(snap_gas, "Distribution spatiale — phase gazeuse", "distrib_gaz.png")
    snap_trans !== nothing && plot_spatial_dist(snap_trans, "Distribution spatiale — transition", "distrib_transition.png")
    snap_solid !== nothing && plot_spatial_dist(snap_solid, "Distribution spatiale — phase solide", "distrib_solide.png")
end


# Q14.9 : compare T=10K (solide) et T=40K (gaz/liquide) dans le même dossier
function neon_lj()
    export_path = "export/q14"
    isdir(export_path) && rm(export_path, recursive=true)
    mkpath(export_path)
    neon_lj_sim(10.0, export_path)
    neon_lj_sim(40.0, export_path)
end


function launch_simulation()
    simple_collision()      # Q4.4
    single_collision()      # Q4.2
    he()                    # Q7.1–Q7.5
    gravity_time()          # Q8.3
    he_with_g()             # Q8.4, Q8.6
    he_ar()                 # Q9.1, Q9.2, Q9.4
    stability()             # Q10.2
    shannon_open_wall()     # Q12.2 + Q12.3
    temp_gradient_sim()     # Q13.3, Q13.4, Q13.5
    neon_lj()               # Q14.7–Q14.9
    condensation_sim()      # Q15.2, Q15.3, Q15.4
end
