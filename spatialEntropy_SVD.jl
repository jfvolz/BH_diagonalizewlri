# Calculate the spatial entanglement entropy of a region A, using the SVD.
function SpatialEE_SVD(N, M, Asize, d)
	# N, M, size of region A, vector of eigenvalues from Lanczos

	Bsize = M - Asize

	# All you need is the sizes of the reduced Hilbert spaces
	DimA = 1 # start at 1, for the [0,0,0,....,0] element
	DimB = 1
	for i=0:N-1
		DimA += binomial(N-i+Asize-1, Asize-1)
		DimB += binomial(N-i+Bsize-1, Bsize-1)
	end

	# dimension of the total Hilbert space
	D = size(basis, 2)

	Amatrix = zeros(DimA, DimB) # This is the matrix that we will SVD

	# form the Amatrix
	for i=1:D
		bra = basis[:, i]
		braA = sub(bra, 1:Asize)
		braB = sub(bra, Asize+1:M)

		row = subSerialNum(sum(braA), Asize, N, braA)
		col = subSerialNum(sum(braB), Bsize, N, braB)

		Amatrix[row, col] = d[2][i] # Assign the matrix the appropriate element from PSI
	end

	U, S, V = svd(Amatrix) # singular value decomposition

	# sum squares of singular values
	sum2=0.0
	sum4=0.0
	for i=1:size(S,1)
		sum2 += S[i]^2
		sum4 += S[i]^4
	end
	if abs(sum2 - 1.0) > 1e-12 # some bounds on the eigenvalue sum
		warn("RDM eigenvalue error ", sum2)
	end

	-log(sum4)
end
