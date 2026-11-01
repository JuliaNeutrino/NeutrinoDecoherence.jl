"""
Environmental decoherence of neutrino oscillations for GLoBES.jl: damping
of the coherences between propagation eigenstates, as a propagation block
that combines with the other blocks in `NewPhysics`.
"""
module NeutrinoDecoherence

using GLoBES
using GLoBES: PropagationBlock
import GLoBES: block_parameters, parameter_aliases, damps, damping_rates

export Decoherence

"""
    Decoherence(; n=0, E0=1.0, pairs=false)

Environmental decoherence: the coherence between the propagation eigenstates
i, j (of the Hamiltonian in matter, numbered by increasing eigenvalue)
decays as exp(-Γ_ij (E/E₀)ⁿ L), with Γ in GeV and E and E₀ in GeV. One
parameter `Γ` for all pairs, or `Γ21`, `Γ31` and `Γ32` with `pairs = true`
(ASCII `Gamma`, `Gamma21`, ...).

The power n sets the energy dependence; n = 0, 2 and -1 are the usual
choices (Lisi, Marrone and Montanino, hep-ph/0002053).

```julia
Decoherence()   # one Γ for all pairs, independent of the energy
Decoherence(n = 2, pairs = true)   # Γ21, Γ31 and Γ32, times (E/E0)² with E0 = 1 GeV
```
"""
struct Decoherence <: PropagationBlock
    n::Float64
    E0::Float64
    pairs::Bool
    function Decoherence(; n::Real = 0.0, E0::Real = 1.0, pairs::Bool = false)
        E0 > 0 || throw(ArgumentError("E0 must be positive, not $E0"))
        return new(n, E0, pairs)
    end
end

block_parameters(b::Decoherence) = b.pairs ? [:Γ21, :Γ31, :Γ32] : [:Γ]
parameter_aliases(b::Decoherence) =
    b.pairs ? [:Gamma21 => :Γ21, :Gamma31 => :Γ31, :Gamma32 => :Γ32] : [:Gamma => :Γ]
damps(::Decoherence) = true
function damping_rates(b::Decoherence, p, x, E)
    f = (E / b.E0)^b.n
    return b.pairs ? (f * x[:Γ21], f * x[:Γ31], f * x[:Γ32]) : ntuple(_ -> f * x[:Γ], 3)
end

Base.show(io::IO, b::Decoherence) = print(
    io,
    "Decoherence(n = $(b.n)",
    b.E0 == 1 ? "" : ", E0 = $(b.E0)",
    b.pairs ? ", pairs = true" : "",
    ")",
)

end
