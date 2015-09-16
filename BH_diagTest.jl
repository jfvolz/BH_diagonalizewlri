#a JULIA program for a full diagonalizaiton 
# for an OBC chain in 1D

#Number of particles
const N = 4 
#Number of sites
const M = 4 

#Hamiltonian parameters
const T = -1.0
const U = 1.0

reload("BH_sparseHam.jl")
#http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
d = eigs(SparseHam, nev=1, which=:SR) 
println(U," ",d[2])

#full diag
#reload("BH_fullHam.jl")
#d = eig(FullHam)
#println(d[2][:,1])
