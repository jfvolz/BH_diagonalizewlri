# Julia program for a full diagonalization for an OBC chain in 1D.

# Number of particles
const N = 4
# Number of sites
const M = 4
# size of region A
const Asize = 2

include("BH_basis.jl")
include("BH_sparseHam.jl")
include("particleEntropy_SVD.jl")
include("spatialEntropy_SVD.jl")

basis = CreateBasis(N, M)

# Hamiltonian parameters
T = -1.0

open("output.dat", "w") do f
	for U=1.0:0.5:20.0
		# Create the Hamiltonian
		SparseHam = CreateSparseHam(basis, T, U)

		# Perform the Lanczos diagonalization to obtain the lowest eigenvector
		# http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
		d = eigs(SparseHam, nev=1, which=:SR)

		# Calculate the second Renyi entropy
		s2_particle = ParticleEE_SVD(N, M, Asize, d[2])
		s2_spatial = SpatialEE_SVD(N, M, Asize, d[2])

		write(f, join((U, d[1][1], s2_particle, s2_spatial), " "), "\n")
		flush(f)
	end
end
