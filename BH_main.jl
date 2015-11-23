# Renyi entanglement entropy of Bose-Hubbard chains in 1D.

using ArgParse

s = ArgParseSettings()
s.autofix_names = true
@add_arg_table s begin
	"M"
		help = "number of sites"
		arg_type = Int
		required = true
	"N"
		help = "number of particles"
		arg_type = Int
		required = true
	"--out"
		metavar = "FILE"
		help = "path to output file"
		required = true
end
add_arg_group(s, "boundary conditions")
@add_arg_table s begin
	"--pbc"
		help = "periodic boundary conditions (default)"
		action = :store_const
		dest_name = "boundary"
		constant = :PBC
	"--obc"
		help = "open boundary conditions"
		action = :store_const
		dest_name = "boundary"
		constant = :OBC
end
add_arg_group(s, "entanglement entropy")
@add_arg_table s begin
	"--ee-all"
		metavar = "XA"
		help = "compute all EEs"
		arg_type = Int
		required = true
end
c = parsed_args = parse_args(ARGS, s, as_symbols=true)

# Number of sites
const M = c[:M]
# Number of particles
const N = c[:N]
# Output file
const output = c[:out]
# Boundary conditions
const boundary = c[:boundary] === nothing ? :PBC : c[:boundary]
# Size of region A
const Asize = c[:ee_all]

include("BH_basis.jl")
include("BH_sparseHam.jl")
include("particleEntropy_SVD.jl")
include("spatialEntropy_SVD.jl")
include("operationalEntropy_SVD.jl")

basis = CreateBasis(N, M)

# Hamiltonian parameters
T = -1.0

open(output, "w") do f
	write(f, "# M=$(M), N=$(N), $(boundary)\n")
	write(f, "# U/t E0 S2(n=$(Asize)) S2(l=$(Asize)) Eop(l=$(Asize))\n")

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
