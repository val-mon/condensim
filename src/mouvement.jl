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
            label=mol.formula,
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
    gif(anim, "export/mouvement.gif", fps=30)
end

function simulate_mouvement()
    single_mol_sim()
end
