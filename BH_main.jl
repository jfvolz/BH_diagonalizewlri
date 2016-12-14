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
    "--site-max"
        metavar = "N"
        help = "site occupation restriction"
        arg_type = Int
end
add_arg_group(s, "output settings")
@add_arg_table s begin
    "--out"
        metavar = "FILE"
        help = "path to output file"
        required = true
    "--no-progress"
        help = "hide progress bar"
        action = :store_true
end
add_arg_group(s, "boundary conditions")
@add_arg_table s begin
    "--pbc"
        help = "periodic boundary conditions (default)"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = PBC
        default = PBC
    "--obc"
        help = "open boundary conditions"
        arg_type = BdryCond
        action = :store_const
        dest_name = "boundary"
        constant = OBC
        default = PBC
end
add_arg_group(s, "BH parameters")
@add_arg_table s begin
    "--u-min"
        metavar = "U"
        help = "minimum U"
        arg_type = Float64
        default = 1.0
    "--u-max"
        metavar = "U"
        help = "maximum U"
        arg_type = Float64
        default = 20.0
    "--u-step"
        metavar = "U"
        help = "U step"
        arg_type = Float64
    "--u-num"
        metavar = "N"
        help = "number of U"
        arg_type = Int
    "--u-log"
        help = "use logarithmic scale for U"
        action = :store_true
    "--t"
        metavar = "t"
        help = "t value"
        arg_type = Float64
        default = 1.0
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

if c[:no_progress]
    include("utils/fakeprogress.jl")
else
    try
        using ProgressMeter
    catch
        include("utils/fakeprogress.jl")
        warn("Install ProgressMeter for a progress bar or use --no-progress to silence this message.")
    end
end

# Number of sites
const M = c[:M]
# Number of particles
const N = c[:N]
# Output file
const output = c[:out]
# Site occupation restriction
const site_max = c[:site_max]
# Boundary conditions
const boundary = c[:boundary]
# Size of region A
const Asize = c[:ee]

if c[:u_log] && c[:u_num] === nothing
    println("--u-log must be used with --u-num")
    exit(1)
end

if c[:u_step] === nothing
    if c[:u_num] === nothing
        U_range = c[:u_min]:0.5:c[:u_max]
    else
        if c[:u_log]
            U_range = logspace(c[:u_min], c[:u_max], c[:u_num])
        else
            U_range = linspace(c[:u_min], c[:u_max], c[:u_num])
        end
    end
else
    if c[:u_num] === nothing
        U_range = c[:u_min]:c[:u_step]:c[:u_max]
    else
        println("--u-step and --u-num may not both be supplied")
        exit(1)
    end
end

if site_max === nothing
    const basis = Szbasis(M, N)
else
    const basis = RestrictedSzbasis(M, N, site_max)
end

open(output, "w") do f
    if site_max === nothing
        write(f, "# M=$(M), N=$(N), $(boundary)\n")
    else
        write(f, "# M=$(M), N=$(N), max=$(site_max), $(boundary)\n")
    end
    write(f, "# U/t E0/t S2(n=$(Asize)) S2(l=$(Asize)) Eop(l=$(Asize))\n")

    @showprogress for U in U_range
        # Create the Hamiltonian
        H = sparse_hamiltonian(basis, c[:t], U, boundary=boundary)

        # Perform the Lanczos diagonalization to obtain the lowest eigenvector
        # http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
        d = eigs(H, nev=1, which=:SR)
        wf = vec(d[2])

        # Calculate the second Renyi entropy
        s2_particle = particle_entropy(basis, Asize, wf)
        s2_spatial, s2_operational = spatial_entropy(basis, Asize, wf)

        write(f, "$(U/c[:t]) $(d[1][1]/c[:t]) $(s2_particle) $(s2_spatial) $(s2_operational)\n")
        flush(f)
    end
end
