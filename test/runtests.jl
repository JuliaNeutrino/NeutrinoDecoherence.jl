using NeutrinoDecoherence, GLoBES, Test

const DUNE = (; L = 1300, ρ = 2.848)
const E = [0.5, 1.0, 2.5, 5.0]

@testset "parameters" begin
    @test NewPhysics(Decoherence()).names == [:Γ]
    m = NewPhysics(Decoherence(n = 2, pairs = true))
    @test m.names == [:Γ21, :Γ31, :Γ32]
    q = OscParams(m, nufit_params("6.1"); Gamma21 = 1e-23)
    @test (q.Γ21, q.Γ31) == (1e-23, 0.0)
    @test_throws ArgumentError Decoherence(E0 = 0)
    @test sprint(show, Decoherence()) == "Decoherence(n = 0.0)"
    @test sprint(show, Decoherence(n = 2, E0 = 10, pairs = true)) ==
          "Decoherence(n = 2.0, E0 = 10.0, pairs = true)"
end

@testset "probabilities" begin
    p = nufit_params("6.1")
    standard = [probability(p, :μ, :e; E, DUNE..., antineutrino = a) for a in (false, true)]
    m = NewPhysics(Decoherence())
    use_model!(m)

    # Γ = 0: the standard probabilities
    @test probability(OscParams(m, p), :μ, :e; E, DUNE...) ≈ standard[1] rtol = 1e-12
    @test probability(OscParams(m, p), :μ, :e; E, DUNE..., antineutrino = true) ≈ standard[2] rtol =
        1e-12

    # Large Γ damps the oscillations to the average (in vacuum, the mass basis)
    U = pmns_matrix(p)
    average = sum(abs2(U[2, i]) * abs2(U[1, i]) for i = 1:3)
    @test probability(OscParams(m, p; Γ = 1e-18), :μ, :e; E = [1.0], L = 1300, ρ = 0.0)[1] ≈ average rtol =
        1e-6

    # Probability is conserved: decoherence only removes interference
    q = OscParams(m, p; Γ = 3e-23)
    total = sum(probability(q, :μ, β; E, DUNE...) for β in (:e, :μ, :τ))
    @test total ≈ ones(length(E)) rtol = 1e-12

    # Equal Γ_ij per pair are one Γ
    per_pair = NewPhysics(Decoherence(pairs = true))
    one_Γ = probability(q, :μ, :e; E, DUNE...)
    use_model!(per_pair)
    @test probability(
        OscParams(per_pair, p; Γ21 = 3e-23, Γ31 = 3e-23, Γ32 = 3e-23),
        :μ,
        :e;
        E,
        DUNE...,
    ) ≈ one_Γ rtol = 1e-12

    # Energy dependence: Γ (E/E0)ⁿ with Γ = Γ0 E0ⁿ at n = 0
    powered = NewPhysics(Decoherence(n = 2, E0 = 2.0))
    use_model!(powered)
    for Eν in E
        P = probability(OscParams(powered, p; Γ = 3e-23), :μ, :e; E = [Eν], DUNE...)[1]
        use_model!(m)
        @test P ≈ probability(OscParams(m, p; Γ = 3e-23 * (Eν / 2)^2), :μ, :e; E = [Eν], DUNE...)[1] rtol =
            1e-12
        use_model!(powered)
    end
    standard_engine!()
end

# P(μ → e), P(μ → μ) of the C engine of SNSI-Decoherence-DUNE (scalar NSI and
# decoherence together): config, sign of CP, E, P_μe, P_μμ
using NeutrinoScalarNSI

const C_REFERENCE = [
    parse.(Float64, split(l)) for
    l in eachline(joinpath(@__DIR__, "snsi_c_reference.txt")) if !startswith(l, "#")
]

c_base(config) =
    config < 3 ?
    OscParams(
        θ12 = asin(sqrt(0.3088)),
        θ13 = asin(sqrt(0.02249)),
        θ23 = asin(sqrt(0.470)),
        δ = deg2rad(207),
        Δm21 = 7.537e-5,
        Δm31 = 2.521e-3,
    ) :
    OscParams(
        θ12 = asin(sqrt(0.3088)),
        θ13 = asin(sqrt(0.02261)),
        θ23 = asin(sqrt(0.555)),
        δ = deg2rad(283),
        Δm21 = 7.537e-5,
        Δm31 = -2.5e-3 + 7.537e-5,
    )

@testset "C engine, with scalar NSI" begin
    m = NewPhysics(ScalarNSI(:eμ, :eτ), Decoherence())
    use_model!(m)
    couplings = Dict(
        1 => (η_eτ = 0.1, Γ = 1e-23),
        2 => (η_eμ = 0.05, δ_eμ = 1.0),
        3 => (η_eτ = -0.05, δ_eτ = π / 3, Γ = 5e-24),
    )
    for (config, cp, Eν, pμe, pμμ) in C_REFERENCE
        q = OscParams(m, c_base(Int(config)); couplings[Int(config)]...)
        anti = cp < 0
        @test probability(q, :μ, :e; E = [Eν], DUNE..., antineutrino = anti)[1] ≈ pμe rtol = 1e-10
        @test probability(q, :μ, :μ; E = [Eν], DUNE..., antineutrino = anti)[1] ≈ pμμ rtol = 1e-10
    end
    standard_engine!()
end
