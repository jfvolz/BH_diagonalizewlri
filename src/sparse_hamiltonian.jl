"""
Number of links for the boundary conditions.
"""
num_links(basis::AbstractSzbasis, boundary::BdryCond) = boundary == PBC ? basis.K : (basis.K - 1)

"""
Create a sparse Hamiltonian matrix for a PBC/OBC BH chain in 1D.

    H = -\\sum_{<i, j>} t_{i,j} (b_i^\\dagger b_j + b_i b_j^\\dagger) + (U/2) \\sum_i n_i (n_i - 1) - \\sum_i \\mu_i n_i
"""
# Function that defines long range potential

function V(r::Int64, Vexponent::Int64, Vscale::Int64)
    Vscale/(r^Vexponent)
end


function sparse_hamiltonian(basis::AbstractSzbasis, Ts::AbstractVector{Float64}, mus::AbstractVector{Float64}, U::Float64, Vexponent::Int64, Vscale::Int64; boundary::BdryCond=PBC)
    end_site = num_links(basis, boundary)

    length(Ts) == end_site || error("Incorrect number of Ts: $(length(Ts)) != $(end_site)")
    length(mus) == basis.K || error("Incorrect number of mus: $(length(mus)) != $(basis.K)")

    rows = Int64[]
    cols = Int64[]
    elements = Float64[]

    for (i, bra) in enumerate(basis)
        # Diagonal part
        Usum = 0
        musum = 0.0
        for j in 1:basis.K
            for r in 0:(basis.K-j)
                if r == 0
                    Usum += U / 2.0 * bra[j] * (bra[j]-1)
                else
                    Usum += V(r, Vexponent, Vscale) * bra[j] * bra[j+r]
                end
            musum += mus[j] * bra[j]
          end
        end
        push!(rows, i)
        push!(cols, i)
        push!(elements, Usum - musum)

        # Off-diagonal part
        for j in 1:end_site
            j_next = j%basis.K + 1
            # Tunnel right, tunnel left.
            for (site1, site2) in [(j, j_next), (j_next, j)]
                if bra[site1] > 0
                    ket = copy(bra)
                    ket[site1] -= 1
                    ket[site2] += 1
                    if ket in basis
                        push!(rows, i)
                        push!(cols, serial_num(basis, ket))
                        push!(elements, -Ts[j] * sqrt(bra[site1]) * sqrt(bra[site2]+1))
                    end
                end
            end
        end
    end

    sparse(rows, cols, elements, length(basis), length(basis))
end



function sparse_hamiltonian(basis::AbstractSzbasis, Ts::AbstractVector{Float64}, U::Float64, Vexponent::Int64, Vscale::Int64; boundary::BdryCond=PBC)
    sparse_hamiltonian(basis, Ts, zeros(basis.K), U, Vexponent, Vscale, boundary=boundary)
end

function sparse_hamiltonian(basis::AbstractSzbasis, T::Float64, mus::AbstractVector{Float64}, U::Float64, Vexponent::Int64, Vscale::Int64; boundary::BdryCond=PBC)
    sparse_hamiltonian(basis, fill(T, num_links(basis, boundary)), mus, U, Vexponent, Vscale, boundary=boundary)
end

function sparse_hamiltonian(basis::AbstractSzbasis, T::Float64, U::Float64, Vexponent::Int64, Vscale::Int64; boundary::BdryCond=PBC)
    sparse_hamiltonian(basis, fill(T, num_links(basis, boundary)), U, Vexponent, Vscale, boundary=boundary)
end
