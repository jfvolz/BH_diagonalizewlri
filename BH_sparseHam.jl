# Create the full Hamiltonian matrix for a PBC/OBC chain in 1D.
# H = -T \sum_<ij> (b+ib-j + b-ib+j) + U/2 \sum_i n_i(n_i-1)
function CreateSparseHam(basis, T, U; boundary=:PBC)
	# Boundary conditions should be :PBC or :OBC.

	I = Int64[] # empty arrays for sparse Hamiltonian
	J = Int64[]
	Element = Float64[]

	for (i, bra) in enumerate(basis)
		# Diagonal part
		Usum = 0
		for j=1:M
			Usum += bra[j] * (bra[j]-1)
		end
		# Sparse Matrix creation
		push!(I, i) # row
		push!(J, i) # column
		push!(Element, U * Usum/2.) # Diagonal operators

		# Off-diagonal part
		end_site = boundary == :PBC ? M : M - 1

		for j=1:end_site
			site1 = j
			if j != M
				site2 = j+1
			else
				site2 = 1
			end

			ket1 = copy(bra)
			A1 = ket1[site1]-1 # We're going to check for annihilation of the state
			A2 = ket1[site2]-1

			# A^dagger A
			if A2 >=0
				ket1[site1] += 1
				ket1[site2] -= 1
				if ket1 in basis
					val1 = sqrt(bra[site1]+1) * sqrt(bra[site2]) # sqrt of occupation
					# Now find the position of the kets using their Serial Number
					b = serial_num(basis, ket1)
					push!(I, i) # row
					push!(J, b) # column
					push!(Element, T * val1)
				end
			end

			ket2 = copy(bra)
			# A A^dagger
			if A1 >=0
				ket2[site1] -= 1
				ket2[site2] += 1
				if ket2 in basis
					val2 = sqrt(bra[site1]) * sqrt(bra[site2]+1) # sqrt of occupation
					# Now find the position of the kets using their Serial Number
					b = serial_num(basis, ket2)
					push!(I, i) # row
					push!(J, b) # column
					push!(Element, T * val2)
				end
			end
		end
	end

	sparse(I, J, Element, length(basis), length(basis)) # create the actual sparse matrix
end
