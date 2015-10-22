# Compute the serial number which indicates the location of a basis vector (see
# equation 12 of http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf).
function SerialNum(N, M, vec) # vec should have length M: i.e. the basis element
	II = 1

	for mu=1:M
		Smu = 0
		for nu=mu+1:M
			Smu += vec[nu]
		end

		for i=0:vec[mu]-1
			II += binomial(N-Smu-i+mu-2, mu-2) # typo: sum missing in paper
		end
	end

	II
end


# Compute the serial number for a partially traced basis.
function subSerialNum(N, M, Ntot, vec)
	offset = 0

	for i=1:Ntot-N # add the dims of Hilbert spaces of other particle numbers
		offset += binomial(Ntot-i+M, M-1)
	end

	SerialNum(N, M, vec) + offset
end
