using Test
using LinearAlgebra
using Condensim

const ATOL = 1e-10
momentum(mol1, mol2) = mol1.mass .* mol1.velocity .+ mol2.mass .* mol2.velocity
energy(mol1, mol2) = 0.5 * mol1.mass * norm(mol1.velocity)^2 + 0.5 * mol2.mass * norm(mol2.velocity)^2
speed(mol) = norm(mol.velocity)

@testset "collisions" begin
    @testset "frontal shock, same masses" begin
        mol1 = Molecule(1.0, 0.1, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "A")
        mol2 = Molecule(1.0, 0.1, [0.2, 0.0, 0.0], [-1.0, 0.0, 0.0], "B")
        p0 = momentum(mol1, mol2)

        resolve_collision!(mol1, mol2)

        @test isapprox(mol1.velocity, [-1.0, 0.0, 0.0]; atol=ATOL)
        @test isapprox(mol2.velocity, [1.0, 0.0, 0.0]; atol=ATOL)
        @test isapprox(momentum(mol1, mol2), p0; atol=ATOL)
    end

    @testset "frontal shock, diff masses" begin
        mol1 = Molecule(4.0, 0.1, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "A")
        mol2 = Molecule(20.0, 0.1, [0.2, 0.0, 0.0], [0.0, 0.0, 0.0], "B")
        p0 = momentum(mol1, mol2)
        e0 = energy(mol1, mol2)

        resolve_collision!(mol1, mol2)

        @test isapprox(mol1.velocity, [-2 / 3, 0.0, 0.0]; atol=ATOL)
        @test isapprox(mol2.velocity, [1 / 3, 0.0, 0.0]; atol=ATOL)
        @test isapprox(momentum(mol1, mol2), p0; atol=ATOL)
        @test isapprox(energy(mol1, mol2), e0; atol=ATOL)
    end

    @testset "45 degree shock" begin
        s2 = sqrt(2) / 2
        mol1 = Molecule(1.0, 0.1, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "A")
        mol2 = Molecule(1.0, 0.1, [0.2 * s2, 0.2 * s2, 0.0], [0.0, 0.0, 0.0], "B")
        p0 = momentum(mol1, mol2)
        e0 = energy(mol1, mol2)

        resolve_collision!(mol1, mol2)

        @test isapprox(mol1.velocity, [0.5, -0.5, 0.0]; atol=ATOL)
        @test isapprox(mol2.velocity, [0.5, 0.5, 0.0]; atol=ATOL)
        @test isapprox(momentum(mol1, mol2), p0; atol=ATOL)
        @test isapprox(energy(mol1, mol2), e0; atol=ATOL)
    end

    @testset "no contact" begin
        mol1 = Molecule(1.0, 0.1, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], "A")
        mol2 = Molecule(1.0, 0.1, [2.0, 0.0, 0.0], [0.0, 0.0, 0.0], "B")
        @test !are_colliding(mol1, mol2)
    end
end

@testset "wall reflections" begin
    # Domain(2,2,2) → walls at ±1, effective boundary at ±(1 - radius) = ±0.9

    @testset "right x wall" begin
        mol = Molecule(1.0, 0.1, [0.95, 0.0, 0.0], [1.0, 0.0, 0.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        s0 = speed(mol)
        reflect_walls!(mol, d)
        @test mol.pos[1] < 0.9 + ATOL
        @test isapprox(mol.velocity, [-1.0, 0.0, 0.0]; atol=ATOL)
        @test isapprox(speed(mol), s0; atol=ATOL)
    end

    @testset "left x wall" begin
        mol = Molecule(1.0, 0.1, [-0.95, 0.0, 0.0], [-1.0, 0.0, 0.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        s0 = speed(mol)
        reflect_walls!(mol, d)
        @test mol.pos[1] > -0.9 - ATOL
        @test isapprox(mol.velocity, [1.0, 0.0, 0.0]; atol=ATOL)
        @test isapprox(speed(mol), s0; atol=ATOL)
    end

    @testset "right y wall" begin
        mol = Molecule(1.0, 0.1, [0.0, 0.95, 0.0], [0.0, 1.0, 0.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        s0 = speed(mol)
        reflect_walls!(mol, d)
        @test mol.pos[2] < 0.9 + ATOL
        @test isapprox(mol.velocity, [0.0, -1.0, 0.0]; atol=ATOL)
        @test isapprox(speed(mol), s0; atol=ATOL)
    end

    @testset "left y wall" begin
        mol = Molecule(1.0, 0.1, [0.0, -0.95, 0.0], [0.0, -1.0, 0.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        s0 = speed(mol)
        reflect_walls!(mol, d)
        @test mol.pos[2] > -0.9 - ATOL
        @test isapprox(mol.velocity, [0.0, 1.0, 0.0]; atol=ATOL)
        @test isapprox(speed(mol), s0; atol=ATOL)
    end

    @testset "right z wall" begin
        mol = Molecule(1.0, 0.1, [0.0, 0.0, 0.95], [0.0, 0.0, 1.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        s0 = speed(mol)
        reflect_walls!(mol, d)
        @test mol.pos[3] < 0.9 + ATOL
        @test isapprox(mol.velocity, [0.0, 0.0, -1.0]; atol=ATOL)
        @test isapprox(speed(mol), s0; atol=ATOL)
    end

    @testset "left z wall" begin
        mol = Molecule(1.0, 0.1, [0.0, 0.0, -0.95], [0.0, 0.0, -1.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        s0 = speed(mol)
        reflect_walls!(mol, d)
        @test mol.pos[3] > -0.9 - ATOL
        @test isapprox(mol.velocity, [0.0, 0.0, 1.0]; atol=ATOL)
        @test isapprox(speed(mol), s0; atol=ATOL)
    end

    @testset "corner x and y" begin
        mol = Molecule(1.0, 0.1, [0.95, 0.95, 0.0], [1.0, 1.0, 0.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        s0 = speed(mol)
        reflect_walls!(mol, d)
        @test mol.pos[1] < 0.9 + ATOL
        @test mol.pos[2] < 0.9 + ATOL
        @test isapprox(mol.velocity, [-1.0, -1.0, 0.0]; atol=ATOL)
        @test isapprox(speed(mol), s0; atol=ATOL)
    end

    @testset "no reflection inside" begin
        mol = Molecule(1.0, 0.1, [0.0, 0.0, 0.0], [1.0, 2.0, 3.0], "A")
        d = Domain(2.0, 2.0, 2.0)
        pos0 = copy(mol.pos)
        vel0 = copy(mol.velocity)
        reflect_walls!(mol, d)
        @test mol.pos == pos0
        @test mol.velocity == vel0
    end
end
