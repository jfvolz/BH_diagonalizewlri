#a JULIA program to create the hilbert space of the Bose Hubbard model
# H = -J \sum_<ij> (b+ib-j + b-ib+j) + U/2 \sum_i n_i(n_i-1)

#Number of particles
N = 4 
#Number of sites
M = 4 

#Hamiltonian parameters
J = -1.0
U = 1.0

#Hilbert Space Dimension
D = binomial(N+M-1,M-1)

println(D)

#start populating the basis vector
starter = zeros(Int64,M)
starter[1] = N

v = starter
println(v)

for i = 2:D 
	min = 0
	j = 1
	while min == 0
		if v[j] != 0 
			min = j
		end
		j += 1
	end
	println("min ",min)


	if min == 1
		v[1] -= 1
		v[2] += 1
	elseif min == M
		println("end ")
	else
		v[1] = v[min]-1
		v[min] = 0
		v[min+1] += 1
	end


	println(v)


end

#append!(highest,v)

#println(highest)
#println(v)