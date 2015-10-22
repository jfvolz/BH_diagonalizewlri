# Create the full Hilbert space of the Bose-Hubbard model.
function CreateBasis(N, M)
	# Hilbert space dimension
	D = binomial(N+M-1, M-1)

	v = zeros(Int64, M)
	v[1] = N

	# start populating the basis vector
	basis = Array(Int64, M, D)
	basis[:, 1] = v

	for i=2:D
		if v[1] != 0
			v[1] -= 1
			v[2] += 1
		else
			min = findfirst(v)

			v[1] = v[min]-1
			v[min] = 0
			v[min+1] += 1
		end

		basis[:, i] = v
	end

	if any(v[1:M-1] .!= 0) || v[M] != N
		warn("Basis construction error")
	end

	basis
end
