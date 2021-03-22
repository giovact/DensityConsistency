"""
Computes first and second moments of `` \\exp[Jx_{i}x_{j} + h_{i}x_{i} +h_{j}x_{j} -\\frac12 x^t \\cdot S \\cdot x + x^t \\cdot y] ``

Input:
* `F::PairIsing <:Factor`: Factor node ``\\psi_{ij}``
* `y::Vector{Float64}`: Array of cavity fields ``\\in \\mathbb{R}^{ 2 }``
* `S::Vector{Float64}`: Cavity coupling matrix ``\\in \\mathbb{R}^{  2 \\times 2 }``
"""
function moments(F::IsingPair, y::Vector{Float64}, S::Array{Float64,2})
    mi = tanh(F.hi + y[1] + atanh(tanh(F.J - S[1,2]) * tanh(F.hj + y[2])))
    mj = tanh(F.hj + y[2] + atanh(tanh(F.J - S[1,2]) * tanh(F.hi + y[1])))
    cij = tanh(F.J - S[1,2] + atanh(tanh(F.hj + y[2]) * tanh(F.hi + y[1]))) - mi*mj
    return [mi;mj],[1-mi^2 cij; cij 1-mj^2]
end

"""
Computes first and second moments of `` \\exp[-\\psi.E(x) -\\frac12 x^t \\cdot S \\cdot x + x^t \\cdot y] ``

Input:
* `F::EnergyFun <:Factor`: Factor node ``\\psi_{a}``
* `y::Vector{Float64}`: Array of cavity fields ``\\in \\mathbb{R}^{ |\\partial a| }``
* `S::Vector{Float64}`: Cavity coupling matrix ``\\in \\mathbb{R}^{  |\\partial a| \\times |\\partial a| }``

Output: first and second moments under tilted "a" distribution -> `av, cov`
where `av` is a vector of size ``|\\partial a|``, `cov`  is a ``|\\partial a| \\times |\\partial a|`` matrix }
"""
function moments(F::EnergyFun, y::Vector{Float64}, S::Array{Float64,2})
    n = length(F.idx)
    N = 1<<n
    av = zeros(n)
    cov = zeros(n,n)
    z = 0
    X = [ [2*((i >> k)%2)-1 for k in 0:(n-1)] for i=0:(N-1) ]
    Y = [ -F.E(X[i]) + X[i]'*(y - 0.5*S*X[i]) for i=1:N ]
    m = maximum(Y)
    @inbounds for i=1:N
        px = exp(Y[i]-m)
        @assert(px>=0)
        z += px
        av .+= X[i] * px
        cov .+= (X[i]*X[i]') * px
    end
    av ./= z
    cov ./= z
    cov .-= av*av'
    return av, cov
end

pearsonmap(x) = sqrt(1-x)/(atanh(sqrt(1-x))*x)



"""
Update the parameters of the gaussian factor "a" according to closure.

Input:

* `F::EnergyFun <:Factor`: Factor node ``\\psi_{a}``
* `y::Vector{Float64}`: Array of cavity fields of size ``|\\partial a|``
* `S::Vector{Float64}`: Cavity coupling matrix of size ``|\\partial a| \\times |\\partial a|``
* `μta::Vector{Float64}`: tilted magnetizations, size ``|\\partial a|`` -> gets updated as ref.
* `Σta::Matrix{Float64}`: Correlation matrix , size ``|\\partial a| \\times |\\partial a|`` -> gets updated as ref.
* `closure::Symbol`: closure protocol
* `η::Float64`: interpolation parameter (valid on :DC closure) -> ``η=0`` gives BP fixed points
* `epsclamp::Float64` clamp tilted moments in [-1+epsclamp, 1-epsclamp]
* `epsconv::Float64`: precision convergence
* `epsmom::Float64`: error over moments at iteration t

Output:

* `ynew::Vector{Float64}`: fields of gauss. factor a @iter t+1, size ``|\\partial a|``
* `Snew::Vector{Float64}`: coupling matrix of gauss. factor a @iter t+1, ``|\\partial a| \\times |\\partial a|``
* `epsmom::Float64`: error over moments at iteration t+1


"""
function setclosure!(F::Factor, y::Vector{Float64},
                        S::Array{Float64,2},
                        µta::Vector{Float64},
                        Σta::Array{Float64,2},
                        closure::Symbol,
                        η::Float64,
                        λ::Float64,
                        epsclamp::Float64,
                        epsmom::Float64)


  av,cov = moments(F, y, S)
  epsmom = max(epsmom , update!(μta,av,0.0), update!(Σta,cov,0.0))
    if closure == :DC
        var = sqrt.(pearsonmap.(clamp.(diag(cov), epsclamp, 1-epsclamp)))
        cov .*= var * var'
        cov = η*cov+(1-η)*Diagonal(cov)
        icov = inv(cov + λ * I)
        Snew = icov-S
        ynew = icov * av - y
    elseif closure == :DC2
        cov = η*cov + (1-η) * Diagonal(cov)
        icov = inv(cov + λ * I)
        means = diag(cov).*atanh.(clamp.(av, -1+epsclamp,1-epsclamp))
        Snew = icov - S
        ynew = icov * means - y
    elseif closure == :EP
        av, cov = moments(F,y,S)
        icov = inv(cov + λ*I)
        Snew = icov - S;
        ynew = icov * av - y
    else
        throw("Select a valid closure")
    end
    return ynew,Snew,epsmom
end
