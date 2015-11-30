# Calculate the operational entanglement entropy of a region A, using the SVD.
# This is the "entanglement of particles" introduced by Wiseman and Vaccaro in
# 2003.
function OperationalEE_SVD(N, M, Asize, d)
	# Dimension of the total Hilbert space
	D = length(d)
	Bsize = M - Asize

	# Matrices to SVD
	Amatrices = []
	for i=0:N
		DimA = num_vectors(basis, i, Asize)
		DimB = num_vectors(basis, N-i, Bsize)

		push!(Amatrices, zeros(Float64, DimA, DimB))
	end

	norms = zeros(Float64, N+1)

	for (i, bra) in enumerate(basis)
		braA = sub(bra, 1:Asize)
		braB = sub(bra, Asize+1:M)

		row = serial_num(basis, Asize, sum(braA), braA)
		col = serial_num(basis, Bsize, sum(braB), braB)

		Amatrices[1 + sum(braA)][row, col] = d[i]
		norms[1 + sum(braA)] += d[i]^2
	end

	norm_err = abs(sum(norms) - 1.0)

	if norm_err > 1e-12
		warn("norm error ", norm_err)
	end

	Ss = [svdvals(Amatrix / sqrt(n)) for (Amatrix, n) in zip(Amatrices, norms)]
	errs = [abs(sum(S.^2) - 1.0) for S in Ss]

	if any(errs .> 1e-12)
		warn("RDM eigenvalue error ", maximum(errs))
	end

	S2s = [-log(sum(S.^4)) for S in Ss]
	dot(norms, S2s)
end
