module BoseHubbardDiagonalize

using JeszenszkiBasis

export
    sparse_hamiltonian,
    particle_entropy,
    spatial_entropy

include("sparse_hamiltonian.jl")
include("particle_entropy.jl")
include("spatial_entropy.jl")

end
