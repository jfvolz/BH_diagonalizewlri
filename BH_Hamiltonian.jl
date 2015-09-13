#a JULIA program to create the Hamiltonian matrix 

#Create the basis
reload("BH_basis.jl")

#println(D," ",basis)

FullHam = zeros(D,D)

for i=1:D

	bra = sub(basis,(i-1)*M+1:(i-1)*M+M) #unpack the bra/ket from the total basis vector
	println(bra)

    Usum = 0
	for j=1:M
		Usum += bra[j]*(bra[j]-1)
	end

    FullHam[i,i] = U*Usum/2.  #The diagonal operators

end


println(FullHam)

