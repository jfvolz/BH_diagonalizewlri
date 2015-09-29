#a JULIA program for a full diagonalizaiton 
# for an OBC chain in 1D


#Number of particles
const N = 7 
#Number of sites
const M = 7 
#size of region A
const Asize = 2

include("BH_basis.jl")
include("BH_sparseHam.jl")
include("spatialEntropy.jl")
include("spatialEntropy_SVD.jl")

basis = CreateBasis(N,M)

#Hamiltonian parameters
T = -1.0

f = open("output.dat","w")
#for U=1.0:100.0
U = 1.0
	
	#Create the Hamiltonian
	SparseHam = CreateSparseHam(basis,T,U)

	#http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
	d = eigs(SparseHam, nev=1, which=:SR) 

    #Calcualte the second Renyi entropy
	@time s2 = SpatialEE(N,M,Asize,d)
	@time s2b = SpatialEE_SVD(N,M,Asize,d)
	println(s2," ",s2b)

	write(f,join((U,d[1][1],s2)," "), "\n")
	flush(f)

#end
close(f)

