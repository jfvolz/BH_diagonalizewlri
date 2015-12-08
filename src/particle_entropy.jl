"""
Calculate the particle entanglement entropy for a subset A, using the SVD.
"""
function particle_entropy(basis::AbstractSzbasis, Asize::Int, d::Vector{Float64})
    M = basis.K
    Bsize = basis.N - Asize
    # Dimensions of partition Hilbert spaces
    DimA = M^Asize
    DimB = M^Bsize

    # Matrix to SVD
    Amatrix = zeros(Float64, DimA, DimB)

    fN = factorial(basis.N)
    occupA = Array(Int, M)
    occup = Array(Int, M)

    for i=1:DimA
        fill!(occupA, 0)

        for k=1:Asize
            occupA[1 + div(i - 1, M^(k - 1)) % M] += 1
        end

        for j=1:DimB
            copy!(occup, occupA)

            for k=1:Bsize
                occup[1 + div(j - 1, M^(k - 1)) % M] += 1
            end

            if occup in basis
                norm = 1 / fN
                for x in occup
                    norm *= factorial(x)
                end

                Amatrix[i, j] = sqrt(norm) * d[serial_num(basis, occup)]
            end
        end
    end

    S = svdvals(Amatrix)
    err = abs(sum(S.^2) - 1.0)

    if err > 1e-12
        warn("RDM eigenvalue error ", err)
    end

    -log(sum(S.^4))
end
