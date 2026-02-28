function create_mol()
    # helium
    He = Molecule(
        4.0026 * PHYSICS.u,
        31e-12,
        ones(3),
        ones(3),
        "He"
    )

    # neon
    Ne = Molecule(
        20.1797 * PHYSICS.u,
        38e-12,
        zeros(3),
        zeros(3),
        "Ne"
    )

    # nitrogen
    N2 = Molecule(
        28.0134 * PHYSICS.u,
        65e-12,
        zeros(3),
        zeros(3),
        "N2"
    )

    # oxygen
    O2 = Molecule(
        31.9988 * PHYSICS.u,
        60e-12,
        zeros(3),
        zeros(3),
        "O2"
    )

    return [He, Ne, N2, O2]
end

function compute_next_pos!(mol, dt)
    acc = zeros(3)
    mol.velocity = mol.velocity + acc * dt
    mol.pos = mol.pos + mol.velocity * dt
end

function single_mol_sim()
    mol = create_mol()[1]
    mol.velocity = 0.8 .* (2 .* rand(3) .- 1)

    t = 0.0
    dt = 0.01
    tfinal = 1.0
    bound = max(maximum(abs.(mol.pos)) + maximum(abs.(mol.velocity)) * tfinal, 1.5)
    axis_lims = (-bound, bound)

    anim = @animate while t < tfinal
        scatter3d([mol.pos[1]], [mol.pos[2]], [mol.pos[3]],
            label=false,
            markersize=6,
            title="t = $(round(t, digits=2))",
            xlims=axis_lims,
            ylims=axis_lims,
            zlims=axis_lims,
            aspect_ratio=:equal
        )

        compute_next_pos!(mol, dt)
        t += dt
    end
    gif(anim, "export/single_mol_sim.gif", fps=30)
end

function multiple_mol_sim()
    mol_types = create_mol()
    mol_number = 20

    all_mol = []
    for _ in 1:mol_number
        new_mol = copy(mol_types[rand(1:length(mol_types))]) # choose one of the molecules
        new_mol.velocity = 0.8 .* (2 .* rand(3) .- 1) # generate a velo btw [-0.8, 0.8]
        push!(all_mol, new_mol)
    end

    t = 0.0
    dt = 0.01
    tfinal = 1.0
    bound = maximum(
        maximum(abs.(m.pos)) + maximum(abs.(m.velocity)) * tfinal for m in all_mol
    )
    axis_lims = (-bound, bound)

    anim = @animate while t < tfinal
        p = plot3d(
            title="t = $(round(t, digits=2))",
            xlims=axis_lims,
            ylims=axis_lims,
            zlims=axis_lims,
            aspect_ratio=:equal
        )
        for mol in all_mol
            scatter3d!(p, [mol.pos[1]], [mol.pos[2]], [mol.pos[3]],
                label=false,
                markersize=6,
            )
        end

        for mol in all_mol
            compute_next_pos!(mol, dt)
        end
        t += dt
    end
    gif(anim, "export/multiple_mol_sim.gif")
end

function simulate_mouvement()
    single_mol_sim()
    multiple_mol_sim()
end
