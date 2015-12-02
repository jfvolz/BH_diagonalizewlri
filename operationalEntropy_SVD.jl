# Calculate the operational entanglement entropy of a region A, using the SVD.
# This is the "entanglement of particles" introduced by Wiseman and Vaccaro in
# 2003.
function OperationalEE_SVD(N::Int, M::Int, A, d::Vector{Float64})
	# Dimension of the total Hilbert space
	D = length(d)

	B = setdiff(1:M, A)

	# Matrices to SVD
	Amatrices = []
	for i=0:N
		DimA = num_vectors(basis, i, length(A))
		DimB = num_vectors(basis, N-i, length(B))

		push!(Amatrices, zeros(Float64, DimA, DimB))
	end

	norms = zeros(Float64, N+1)

	for (i, bra) in enumerate(basis)
		braA = sub(bra, A)
		braB = sub(bra, B)

		row = serial_num(basis, length(A), sum(braA), braA)
		col = serial_num(basis, length(B), sum(braB), braB)

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

OperationalEE_SVD(N::Int, M::Int, Asize::Int, d::Vector{Float64}) = OperationalEE_SVD(N, M, 1:Asize, d)
