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

	S = svdvals(Amatrix)
	err = abs(sum(S.^2) - 1.0)

	if err > 1e-12
		warn("RDM eigenvalue error ", err)
	end

	-log(sum(S.^4))
end
