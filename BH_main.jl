# Renyi entanglement entropy of Bose-Hubbard chains in 1D.

push!(LOAD_PATH, joinpath(dirname(@__FILE__), "src"))
using BoseHubbardDiagonalize

using ArgParse
using JeszenszkiBasis

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
    "--site-max"
        metavar = "N"
        help = "site occupation restriction"
        arg_type = Int
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
add_arg_group(s, "BH parameter range")
@add_arg_table s begin
    "--u-min"
        metavar = "U"
        help = "minimum U/t"
        arg_type = Float64
        default = 1.0
    "--u-max"
        metavar = "U"
        help = "maximum U/t"
        arg_type = Float64
        default = 20.0
    "--u-step"
        metavar = "U"
        help = "U/t step"
        arg_type = Float64
        default = 0.5
end
add_arg_group(s, "entanglement entropy")
@add_arg_table s begin
    "--ee"
        metavar = "XA"
        help = "compute all EEs with partition size XA"
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
# Site occupation restriction
const site_max = c[:site_max]
# Boundary conditions
const boundary = c[:boundary] === nothing ? :PBC : c[:boundary]
# Size of region A
const Asize = c[:ee]

if site_max === nothing
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end

# Hamiltonian parameters
T = -1.0

open(output, "w") do f
    if site_max === nothing
        write(f, "# M=$(M), N=$(N), $(boundary)\n")
    else
        write(f, "# M=$(M), N=$(N), max=$(site_max), $(boundary)\n")
    end
    write(f, "# U/t E0 S2(n=$(Asize)) S2(l=$(Asize)) Eop(l=$(Asize))\n")

    for U=c[:u_min]:c[:u_step]:c[:u_max]
        # Create the Hamiltonian
        H = sparse_hamiltonian(basis, T, U, boundary=boundary)

        # Perform the Lanczos diagonalization to obtain the lowest eigenvector
        # http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
        d = eigs(H, nev=1, which=:SR)
        wf = vec(d[2])

        # Calculate the second Renyi entropy
        s2_particle = particle_entropy(basis, Asize, wf)
        s2_spatial, s2_operational = spatial_entropy(basis, Asize, wf)

        write(f, "$(U) $(d[1][1]) $(s2_particle) $(s2_spatial) $(s2_operational)\n")
        flush(f)
    end
end
