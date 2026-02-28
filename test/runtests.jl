using Test
using LinearAlgebra
using Condensim: Molecule, resolve_collision!, are_colliding

const ATOL = 1e-10
momentum(mol1, mol2) = mol1.mass .* mol1.velocity .+ mol2.mass .* mol2.velocity
energy(mol1, mol2) = 0.5 * mol1.mass * norm(mol1.velocity)^2 + 0.5 * mol2.mass * norm(mol2.velocity)^2

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
