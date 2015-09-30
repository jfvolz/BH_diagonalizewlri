# This function calculates the spatial entanglement entropy
# of a region A, using the SVD 
# 

function SpatialEE_SVD(N,M,Asize,d)
#N, M, size of region A, vector of eigenvalues from Lanczos

Bsize = M - Asize

#All you need is the sizes of the reduced hilbert spaces
DimA = 1 #start at 1, for the [0,0,0,....,0] element
DimB = 1
for i=0:N-1
	DimA += binomial(N-i+Asize-1,Asize -1)
	DimB += binomial(N-i+Bsize-1,Bsize -1)
end

#dimension of the total Hilbert space
D = div(length(basis),M)
#println("D ",D)
#dimension of the subspace A
#println("DimA ",DimA)
#println("DimB ",DimB)

Amatrix = zeros(DimA,DimB)  #This is the matrix that we will SVD

#form the Amatrix
for i=1:D

 	bra = sub(basis,(i-1)*M+1:(i-1)*M+M) #unpack the bra/ket from the total basis vector 
	braA = sub(bra,1:Asize)
	braB = sub(bra,Asize+1:M)

	row = subSerialNum(sum(braA),Asize,N,braA)
	col = subSerialNum(sum(braB),Bsize,N,braB)

	Amatrix[row,col] = d[2][i] #Assign the matrix the appropriate element from PSI
end

#println(Amatrix)

(U, S, V) = svd(Amatrix) #singular value decomposition

#println(S,size(S))

#sum squares of singular values
#Ssize = min(DimA,DimB)
#println("Ssize ",Ssize," ",size(S,1))
sum2=0.0
sum4=0.0
for i=1:size(S,1)
	sum2 += S[i]^2
	sum4 += S[i]^4
end
#println(sum2)
if (sum2 - 1.0) > 1e-12  #some bounds on the eigenvalue sum
	println("RDM eigenvalue error ",sum2)
end

#println(-log(sum4))

#return the second Renyi entropy
return -log(sum4)

end

