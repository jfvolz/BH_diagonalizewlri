#a JULIA program to create the full Hamiltonian matrix 
# for a PBC/OBC chain in 1D
# H = -T \sum_<ij> (b+ib-j + b-ib+j) + U/2 \sum_i n_i(n_i-1)

function CreateSparseHam(basis)

#Create the basis
#reload("BH_basis.jl")

I = Int64[]  #empty arrays for Sparse Hamiltonian
J = Int64[] 
Element = Float64[]

D = div(length(basis),M)
println(D)

for i=1:D

	bra = sub(basis,(i-1)*M+1:(i-1)*M+M) #unpack the bra/ket from the total basis vector
	#println(bra)
	  
    #Diagonal part
    Usum = 0
	for j=1:M
		Usum += bra[j]*(bra[j]-1)  
	end
	#Sparse Matrix creation
	push!(I,i) #row
	push!(J,i) #column
	push!(Element,U*Usum/2.)  #The diagonal operators


	#The off-diagonal part 
	#for j=1:M-1 # (OBC chain)
	for j=1:M  # (PBC chain)
		ket1 = copy(bra)
		site1 = j
		if j != M
			site2 = j+1
		else
			site2 = 1
		end

		#Adag1 = ket1[site1]+1
		A1 = ket1[site1]-1
		#Adag2 = ket1[site2]+1
		A2 = ket1[site2]-1

        #A^dagger A
		if A2 >=0
			ket1[site1] += 1
			ket1[site2] -= 1 
			val1 = sqrt(bra[site1]+1) * sqrt(bra[site2]) #sqrt of occupation
		end

		ket2 = copy(bra) 
		#A A^dagger
		if A1 >=0
			ket2[site1] -= 1
			ket2[site2] += 1
			val2 = sqrt(bra[site1]) * sqrt(bra[site2]+1) #sqrt of occupation
		end

		#Now search the entire basis for the corresponding bra
		found1=false
		found2=false
		for b = (i+1):D
			testbra = sub(basis,(b-1)*M+1:(b-1)*M+M) #unpack 
			if testbra == ket1
				push!(I,i) #row
				push!(J,b) #column
				push!(Element,T*val1)  
				push!(I,b) #row
				push!(J,i) #column
				push!(Element,T*val1)  
				found1 = true
			end
			if testbra == ket2
				push!(I,i) #row
				push!(J,b) #column
				push!(Element,T*val2)  
				push!(I,b) #row
				push!(J,i) #column
				push!(Element,T*val2)  
				found2 = true
			end
			if (found1 == true) && (found2 == true) break end
		end
	end
end

SparseHam = sparse(I,J,Element,D,D) #create the actual sparse matrix

return SparseHam
end
