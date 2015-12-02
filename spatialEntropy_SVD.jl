# Calculate the spatial entanglement entropy of a region A, using the SVD.
function SpatialEE_SVD(N::Int, M::Int, A, d::Vector{Float64})
	# N, M, sites in region A, ground state components

	B = setdiff(1:M, A)

	# All you need is the sizes of the reduced Hilbert spaces
	DimA = 1 # start at 1, for the [0,0,0,....,0] element
	DimB = 1
	for i=0:N-1
		DimA += num_vectors(basis, N-i, length(A))
		DimB += num_vectors(basis, N-i, length(B))
	end

	Amatrix = zeros(DimA, DimB) # This is the matrix that we will SVD

	# form the Amatrix
	for (i, bra) in enumerate(basis)
		row = sub_serial_num(basis, sub(bra, A))
		col = sub_serial_num(basis, sub(bra, B))

		Amatrix[row, col] = d[i] # Assign the matrix the appropriate element from PSI
	end

	S = svdvals(Amatrix)
	err = abs(sum(S.^2) - 1.0)

	if err > 1e-12
		warn("RDM eigenvalue error ", err)
	end

	-log(sum(S.^4))
end

SpatialEE_SVD(N::Int, M::Int, Asize::Int, d::Vector{Float64}) = SpatialEE_SVD(N, M, 1:Asize, d)
