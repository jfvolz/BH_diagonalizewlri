# Renyi entanglement entropy of Bose-Hubbard chains in 1D.

if length(ARGS) != 5
    println("usage: <particles> <sites> <A length> <boundary conditions> <output>")
    exit(1)
end

# Number of particles
const N = parse(Int, ARGS[1])
# Number of sites
const M = parse(Int, ARGS[2])
# Size of region A
const Asize = parse(Int, ARGS[3])
# Boundary conditions
const boundary = ARGS[4]
# Output file.
const output = ARGS[5]

if !(boundary in ["OBC", "PBC"])
    println("valid boundary conditions: OBC PBC")
    exit(1)
end

include("BH_basis.jl")
include("BH_sparseHam.jl")
include("particleEntropy_SVD.jl")
include("spatialEntropy_SVD.jl")
include("operationalEntropy_SVD.jl")

basis = CreateBasis(N, M)

# Hamiltonian parameters
T = -1.0

open(output, "w") do f
	for U=1.0:0.5:20.0
		# Create the Hamiltonian
		SparseHam = CreateSparseHam(basis, T, U, boundary=boundary)

		# Perform the Lanczos diagonalization to obtain the lowest eigenvector
		# http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
		d = eigs(SparseHam, nev=1, which=:SR)

		# Calculate the second Renyi entropy
		s2_particle = ParticleEE_SVD(N, M, Asize, d[2])
		s2_spatial = SpatialEE_SVD(N, M, Asize, d[2])
		s2_operational = OperationalEE_SVD(N, M, Asize, d[2])

		write(f, join((U, d[1][1], s2_particle, s2_spatial, s2_operational), " "), "\n")
		flush(f)
	end
end
