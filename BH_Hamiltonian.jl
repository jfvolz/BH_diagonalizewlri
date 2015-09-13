#a JULIA program to create the Hamiltonian matrix 
# for an OBC chain in 1D

#Create the basis
reload("BH_basis.jl")

#println(D," ",basis)

FullHam = zeros(D,D)

for i=1:D

	bra = sub(basis,(i-1)*M+1:(i-1)*M+M) #unpack the bra/ket from the total basis vector
	println(bra)
	  
    #Diagonal part
    Usum = 0
	for j=1:M
		Usum += bra[j]*(bra[j]-1)  
	end
    FullHam[i,i] = U*Usum/2.  #The diagonal operators

	#The off-diagonal part (OBC chain)
	for j=1:M-1
		ket1 = copy(bra)
		site1 = j
		site2 = j+1

		#Adag1 = ket1[site1]+1
		A1 = ket1[site1]-1
		#Adag2 = ket1[site2]+1
		A2 = ket1[site2]-1

        #A^dagger A
		if A2 >=0
			ket1[site1] += 1
			ket1[site2] -= 1
		end

		ket2 = copy(bra) 
		#A A^dagger
		if A1 >=0
			ket2[site1] -= 1
			ket2[site2] += 1
		end

		#Now search the entire basis for the corresponding bra
		found1=false
		found2=false
		for b = (i+1):D
			testbra = sub(basis,(b-1)*M+1:(b-1)*M+M) #unpack 
			if testbra == ket1
				FullHam[i,b] = J
				FullHam[b,i] = J
				found1 = true
			end
			if testbra == ket2
				FullHam[i,b] = J
				FullHam[b,i] = J
				found2 = true
			end
			if found1 == true && found2==true break end
		end
	end
end

println(FullHam)

