#This computes the serial number which indicates the location of a basis vector
# see http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf
# equation 12

function SerialNum(N,M,vec) #vec should have length M: i.e. the basis element

#println("serial#")
#println(vec)
II = 1
for mu=1:M
	#println(mu," ",vec[mu])

	Smu = 0
	for nu = mu+1:M
		#println("nu ",nu)
		Smu += vec[nu]
	end
	#println(mu," Smu2 ",Smu)

    Imu = 0
	for i=0:vec[mu]-1
		Imu += binomial(N-Smu - i + mu - 2,mu-2) #typo: sum missing in paper
	end
	#println(mu," Imu ",Imu)

	II += Imu

end

#println(II)
return II

end

#------------
#This function computes the serial number for a partially traced basis 
function subSerialNum(N,M,Ntot,vec) 

offset = 0
for i = 1:(Ntot-N)        #add the dims of Hilbert spaces of other particle numbers 
	count = Ntot - i + 1  
	offset += binomial(count+M-1,M-1)
end

return SerialNum(N,M,vec) + offset

end



