# NeutrinoDecoherence.jl

Environmental decoherence of neutrino oscillations for [GLoBES.jl](https://github.com/JuliaNeutrino/GLoBES.jl). Interactions with an environment damp the coherence between the propagation eigenstates i, j by exp(-Γ_ij (E/E₀)ⁿ L). This removes interference without changing the total probability.

```julia
using Pkg
Pkg.develop(path = "/path/to/NeutrinoDecoherence.jl")   # not registered yet; Julia 1.11 or later
using GLoBES, NeutrinoDecoherence
```

`Decoherence` is a GLoBES.jl propagation block, with one parameter `Γ` (GeV), or `Γ21`, `Γ31` and `Γ32` with `pairs = true`:

```julia
model = NewPhysics(Decoherence(n = 0))    # Decoherence(n = 2, E0 = 1.0, pairs = true) for Γ_ij (E/E0)²
use_model!(model)                         # the probability engine of all experiments
p = OscParams(model, nufit_params("6.1"); Γ = 1e-23)
probability(p, :μ, :e; E = 2.5, L = 1300, ρ = 2.848)
```

It combines with the other blocks, e.g. `NewPhysics(ScalarNSI(:eτ), Decoherence())` with [NeutrinoScalarNSI.jl](https://github.com/JuliaNeutrino/NeutrinoScalarNSI.jl). Rates, χ² and fits in GLoBES.jl then use Γ like the standard parameters.

Together with scalar NSI, the probabilities agree with the C engine of the SNSI-decoherence DUNE study to 10⁻¹⁰ (`test/snsi_c_reference.txt`). Tests: `julia --project=. test/runtests.jl`.
