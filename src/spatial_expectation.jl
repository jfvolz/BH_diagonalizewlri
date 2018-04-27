"""
calculate <n_i> and <n_i^2> and delta n for a wavefunction
"""
function spatial_expectation{T<:Number}(basis::AbstractSzbasis, wf::Vector{T})
  expectations = zeros(basis.K)
  squareExpectations = zeros(basis.K)
  deltaNsquared = zeros(basis.K)
  for i = 1:basis.K
    nSum = 0
    nSquareSum = 0
    for j = 1:basis.D
      basisState = getindex(basis,j)
      nSum += basisState[i] * abs2(wf[j])
      nSquareSum += basisState[i]^2 * abs2(wf[j])
    end
    expectations[i] += nSum
    squareExpectations[i] += nSquareSum
    deltaNsquared[i] += nSquareSum - nSum^2
  end
  expectations, squareExpectations, deltaNsquared
end
