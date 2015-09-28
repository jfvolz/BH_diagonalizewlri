#a JULIA program for a full diagonalizaiton 
# for an OBC chain in 1D

#Number of particles
const N = 4 
#Number of sites
const M = 4 

#Hamiltonian parameters
const T = -1.0
const U = 1.0

include("BH_basis.jl")
basis = CreateBasis(N,M)
#println(basis)

include("BH_sparseHam.jl")
SparseHam = CreateSparseHam(basis,T,U)
 
#http://docs.julialang.org/en/release-0.3/stdlib/linalg/?highlight=lanczos
d = eigs(SparseHam, nev=1, which=:SR) 
println(U," ",d[2])
  
#------full diag
#reload("BH_fullHam.jl")
#d = eig(FullHam)
#println(d[2][:,1])
#-------

#Begin calculation of the reduced density matrix
Asize = 2 # how many sites in region A
Bsize = M - Asize
Abasis = Int64[] 
Bbasis = Int64[] 
for i=0:N-1
	tempAbasis = CreateBasis(N-i,Asize)
	tempBbasis = CreateBasis(N-i,Bsize)
	append!(Abasis,tempAbasis)
	append!(Bbasis,tempBbasis)
end
tempAbasis = zeros(Asize)
tempBbasis = zeros(Bsize)
append!(Abasis,tempAbasis) #this is the reduced basis for region A
append!(Bbasis,tempBbasis) #this is the reduced basis for region B
println(Abasis)
println(Bbasis)

#dimension of the total Hilbert space
D = div(length(basis),M)
println("D ",D)
#dimension of the subspace A
DimA = div(length(Abasis),Asize)
DimB = div(length(Bbasis),Bsize)
println("DimA ",DimA)
println("DimB ",DimB)

BDensityM = zeros(DimB,DimB) #The density matrix over B
 
#perform the Trace_A
for j=1:DimA
	braA = sub(Abasis,(j-1)*Asize+1:(j-1)*Asize+Asize)	#unpack the bra/ket 
	#index = subSerialNum(sum(braA),Asize,N,braA)		#this is the position in the basis
	#println(j," ",index," ",braA," ",sum(braA))

 	for row=1:D
 		bra = sub(basis,(row-1)*M+1:(row-1)*M+M) #unpack the bra/ket from the total basis vector
 		subbra = sub(bra,1:Asize)
 		#println("trace1 ",braA," ",subbra)
 		#println(bra," ",d[2][row]) #this is the eigenvector 

		for col=1:D
 			ket= sub(basis,(col-1)*M+1:(col-1)*M+M) #unpack the bra/ket from the total basis vector 
			subket = sub(ket,1:Asize)
			#println("trace2 ",braA," ",subket)

			if subbra == braA && subket == braA #orthogonality is satisfied for the trace 
				Bbra = sub(bra,Asize+1:M)
				Bket = sub(ket,Asize+1:M)
				println("Bbra ",braA," ",Bbra," ",Bket)
				indexRow = 	subSerialNum(sum(Bbra),Bsize,N,Bbra)
				indexCol = 	subSerialNum(sum(Bket),Bsize,N,Bket)
				BDensityM[indexRow,indexCol] += d[2][row]*d[2][col]
			end

		end #col
 
	end #row

end #trace

println(BDensityM)
lambda = eig(BDensityM)
println(lambda[1])
println(size(lambda[1]))

Lsum=0.0
for i=1:DimB
	Lsum += lambda[1][i]
end
println(Lsum)

