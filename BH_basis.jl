#a JULIA program to create the full hilbert space of the Bose Hubbard model
# H = -T \sum_<ij> (b+ib-j + b-ib+j) + U/2 \sum_i n_i(n_i-1)


#Hilbert Space Dimension
D = binomial(N+M-1,M-1)
#println(D)

#start populating the basis vector
basis = zeros(Int64,M)
basis[1] = N

v = copy(basis)
#println(v)

i = 1
exit = false
while exit == false   #for i = 2:D 
	i += 1            #Hilbert space counter
	#println("count ",i)

	min = 0
	j = 1
	while min == 0
		if v[j] != 0 
			min = j
		end
		j += 1
	end
	#println("min ",min)

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
		append!(basis,v)
		#println(v)
	end

end

if D != i
	println("Basis dimension error")
end

#println(basis)
