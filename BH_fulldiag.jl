#a JULIA program for a full diagonalizaiton 
# for an OBC chain in 1D

#Create the basis
reload("BH_Hamiltonian.jl")

#http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
d = eigs(FullHam)  
#d = eigvals(FullHam)
println(d)

SparseHam = sparse(FullHam)
println(SparseHam)

d2 = eigs(SparseHam)
println(d2)
