# Create the full Hilbert space of the Bose-Hubbard model.
function CreateBasis(N, M)
	# Hilbert space dimension
	D = binomial(N+M-1, M-1)

	# start populating the basis vector
	basis = zeros(Int64, M)
	basis[1] = N

	v = copy(basis)

	i = 1
	exit = false
	while exit == false # for i = 2:D
		i += 1 # Hilbert space counter

		min = 0
		j = 1
		while min == 0
			if v[j] != 0
				min = j
			end
			j += 1
		end

		if min == M
			i -= 1
			exit = true
		else
			if min == 1
				v[1] -= 1
				v[2] += 1
			else
				v[1] = v[min]-1
				v[min] = 0
				v[min+1] += 1
			end
			append!(basis, v)
		end
	end

	if D != i
		warn("Basis dimension error")
	end

	basis
end
