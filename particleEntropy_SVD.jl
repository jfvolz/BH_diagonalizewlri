# Calculate the particle entanglement entropy for a subset A, using the SVD.
function ParticleEE_SVD(N, M, Asize, d)
	# Dimension of the total Hilbert space
	D = length(d)
	Bsize = N - Asize
	# Dimensions of partition Hilbert spaces
	DimA = M^Asize
	DimB = M^Bsize

	# Matrix to SVD
	Amatrix = zeros(DimA, DimB)

	for i=1:D
		bra = basis[:, i]

		# Generate a representative state from the occupations.
		state = vcat([[j for _=1:bra[j]] for j=1:M]...)
		# Find all the unique permutations of the representative state.
		perms = unique(permutations(state))
		norm = 1.0 / sqrt(length(perms))

		for perm in perms
			stateA = sub(perm, 1:Asize)
			stateB = sub(perm, Asize+1:M)

			row = 1 + sum([M^(i-1) * (x-1) for (i, x) in enumerate(stateA)])
			col = 1 + sum([M^(i-1) * (x-1) for (i, x) in enumerate(stateB)])

			Amatrix[row, col] = norm * d[i]
		end
	end

	S = svdvals(Amatrix)
	err = abs(sum(S.^2) - 1.0)

	if err > 1e-12
		warn("RDM eigenvalue error ", err)
	end

	-log(sum(S.^4))
end
