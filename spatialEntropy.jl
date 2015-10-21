# This function calculates the spatial entanglement entropy
# of a region A
# 

function SpatialEE(N,M,Asize,d)
#N, M, size of region A, vector of eigenvalues from Lanczos

#Begin calculation of the reduced density matrix
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
#println(Abasis)
#println(Bbasis)

#dimension of the total Hilbert space
D = div(length(basis),M)
#println("D ",D)
#dimension of the subspace A
DimA = div(length(Abasis),Asize)
DimB = div(length(Bbasis),Bsize)
#println("DimA ",DimA)
#println("DimB ",DimB)

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
				#println("Bbra ",braA," ",Bbra," ",Bket)
				indexRow = 	subSerialNum(sum(Bbra),Bsize,N,Bbra)
				indexCol = 	subSerialNum(sum(Bket),Bsize,N,Bket)
				BDensityM[indexRow,indexCol] += d[2][row]*d[2][col]
			end

		end #col
 
	end #row

end #trace

#println(BDensityM)
lambda = eig(BDensityM)
#println(lambda[1])
#println(size(lambda[1]))

#sum eigenvalues of RDM
Lsum=0.0
for i=1:DimB
	Lsum += lambda[1][i]
end
if abs(Lsum - 1.0) > 1e-12  #some bounds on the eigenvalue sum
	println("RDM eigenvalue error ",Lsum)
end

#von Neumann entropy
s2 = 0.0
for i=1:DimB
	if lambda[1][i] > 10e-10 
		#println(lambda[1][i])
		#vN -= (lambda[1][i]) * (log(lambda[1][i]) )
		s2 += lambda[1][i] ^2
	end
end
#println("vN entropy ",vN)

return -log(s2)

end

