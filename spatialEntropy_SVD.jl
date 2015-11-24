# Calculate the spatial entanglement entropy of a region A, using the SVD.
function SpatialEE_SVD(N, M, Asize, d)
	# N, M, size of region A, ground state components

	Bsize = M - Asize

	# All you need is the sizes of the reduced Hilbert spaces
	DimA = 1 # start at 1, for the [0,0,0,....,0] element
	DimB = 1
	for i=0:N-1
		DimA += num_vectors(basis, N-i, Asize)
		DimB += num_vectors(basis, N-i, Bsize)
	end

	Amatrix = zeros(DimA, DimB) # This is the matrix that we will SVD

	# form the Amatrix
	for (i, bra) in enumerate(basis)
		braA = sub(bra, 1:Asize)
		braB = sub(bra, Asize+1:M)

		row = sub_serial_num(basis, braA)
		col = sub_serial_num(basis, braB)

		Amatrix[row, col] = d[i] # Assign the matrix the appropriate element from PSI
	end

	S = svdvals(Amatrix)
	err = abs(sum(S.^2) - 1.0)

	if err > 1e-12
		warn("RDM eigenvalue error ", err)
	end

	-log(sum(S.^4))
end
