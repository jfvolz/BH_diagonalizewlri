module BoseHubbardDiagonalize

using JeszenszkiBasis

export
    BdryCond,
    OBC,
    PBC,

    sparse_hamiltonian,
    particle_entropy,
    spatial_entropy

"""
Boundary conditions.
"""
@enum BdryCond PBC OBC
@doc "Periodic boundary conditions." PBC
@doc "Open boundary conditions." OBC

include("sparse_hamiltonian.jl")
include("particle_entropy.jl")
include("spatial_entropy.jl")

end
