"""
Abstract Factor type
"""
abstract type Factor end

"""
Factor type defined by energy function

`` \\psi_{a}(x_{\\partial a}) = exp[ -E(x_{\\partial a})]``

Parameters:
* `E::Function`: ``E(x_{\\partial a}) ``
* `idx::Vector{Int64}`: ``\\partial a`` indeces

"""
struct EnergyFun <: Factor
    E::Function
    idx::Vector{Int64}
end

"""
Factor type for Ising model

`` \\psi_{ij}(x_{i},x_{j}) = exp[ Jx_{i}x_{j} + h_{i}x_{i} +h_{j}x_{j}]``

Parameters:
* `J::Float64`: coupling
* `hi::Float64`: field on ``x_{i}``
* `hi::Float64`: field on ``x_{j}``
* `idx::Vector{Int64}`: ``(i,j)``

moments are computed analytically

"""
struct IsingPair <: Factor
    J::Float64
    hi::Float64
    hj::Float64
    idx::Vector{Int64}
	IsingPair(J,hi,hj,idx) = length(idx)!=2 ? error("idx must have length 2") : new(J,hi,hj,idx)
end

"""
Basic constructor for Ising Factor Graph

Input:
* `N::Int64`: number of spins
* `h::Vector{Float64}`: Array of fields ``\\in \\mathbb{R}^{N}``
* `J::SparseMatrixCSC{Float64,Int64}`: symmetric coupling matrix ``\\in \\mathbb{R}^{N\\times N}``
* `t::Symbol`: choose between `EnergyFun` and `IsingPair`

Output: `Array{Factor}`: array of factor nodes `` \\{ \\psi_{ij} \\}``

"""
function IsingFG(N::Int64,h::Vector{Float64},J::SparseMatrixCSC{Float64,Int64}, t::Symbol; Adj::SparseMatrixCSC{Int64,Int64} = (J.!=0)*1)
	t ∉ (:EnergyFun,:IsingPair) && error("invalid option")
	degree = sum(Adj,dims = 2)[:]
	if t == :EnergyFun
		return collect([EnergyFun((x::Vector{Int64}-> -( J[i,j]*x[1]*x[2]+x[1]*(h[i]/degree[i])+ x[2]*(h[j]/degree[j]))),[i,j] ) for i=1:N, j=1:N if Adj[i,j]!=0 && i<j])
	elseif t == :IsingPair
		return collect([IsingPair(J[i,j],h[i]/degree[i],h[j]/degree[j],[i;j] ) for i=1:N, j=1:N if Adj[i,j]!=0 && i<j])
	end
end

"""
DCState

Istantaneous state of Density Consistency run

Fields:
* `μtl :: Vector{Vector}`		 array of first moments of tilted distributions (one for each factor)
* `Σtl :: Vector{Matrix}`		 array of second moments of tilted distributions (one for each factor)
* `μt  :: Vector`				 mean vector w.r.t tilted distribution (size N)
* `Σt  :: Matrix`				 correlation matrix w.r.t tilted distribution (size NxN)
* `µ   :: Vector`				 mean vector of gaussian distribution (size N)
* `Σ   :: Matrix`				 correlation matrix of the gaussian distribution  (size NxN)
* `h   :: Vector{Vector}`		 array of gaussian field vector (one for each factor)
* `S   :: Vector{Matrix}`		 array of inverse correlation matrix (one for each factor)
* `yc  :: Vector{Vector}`		 array of cavity fields vector (one for each factor)
* `Sc  :: Vector{Matrix}`		 array of cavity inverse correlation matrix (one for each factor)
"""
struct DCState
    μtl :: Vector{Vector}		# array of first moments of tilted distributions (one for each factor)
    Σtl :: Vector{Matrix}		# array of second moments of tilted distributions (one for each factor)
    μt  :: Vector			# mean vector w.r.t tilted distribution (size N)
    Σt  :: Matrix			# correlation matrix w.r.t tilted distribution (size NxN)
    µ   :: Vector			# mean vector of gaussian distribution (size N)
    Σ   :: Matrix			# correlation matrix of the gaussian distribution  (size NxN)
    h   :: Vector{Vector}		# array of gaussian field vector (one for each factor)
    S   :: Vector{Matrix}		# array of inverse correlation matrix (one for each factor)
    yc  :: Vector{Vector}		# array of cavity fields vector (one for each factor)
    Sc  :: Vector{Matrix}		# array of cavity inverse correlation matrix (one for each factor)
end

function DCState(N::Int,m::Vector{Int};h0::Float64 = 1e-3,f::Function=n->h0*(rand(n) .- 0.5))
    M=length(m)
    return DCState([zeros(m[a]) for a=1:M],
                   [zeros(m[a],m[a]) for a=1:M],
                   zeros(N),
                   zeros(N,N),
				   zeros(N),
                   zeros(N,N),
                   [f(m[a]) for a=1:M],
                   [eye(m[a]) for a=1:M],
                   [zeros(m[a]) for a=1:M],
                   [zeros(m[a],m[a]) for a=1:M])
end
