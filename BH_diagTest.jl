#a JULIA program for a full diagonalizaiton 
# for an OBC chain in 1D

#Number of particles
const N = 2 
#Number of sites
const M = 3 

#Hamiltonian parameters
const T = -1.0
const U = 1.0

reload("BH_sparseHam.jl")
#http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
d = eigs(SparseHam, nev=1, which=:SR) 
println(U," ",d[2])

#------full diag
#reload("BH_fullHam.jl")
#d = eig(FullHam)
#println(d[2][:,1])
#-------

const Asize = 1 # how many sites in region A

for i=1:D
	bra = sub(basis,(i-1)*M+1:(i-1)*M+M) #unpack the bra/ket from the total basis vector
	println(bra," ",d[2][i]) #this is the eigenvector 
end

