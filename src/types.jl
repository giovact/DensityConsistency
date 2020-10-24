abstract type Factor end

struct EnergyFun <: Factor
    E::Function
    idx::Vector{Int64}
end

struct PairIsing <: Factor
    J::Float64
    hi::Float64
    hj::Float64
    idx::Vector{Int64}
end

function IsingFactors(N::Int64,h::Vector{Float64},J::SparseMatrixCSC{Float64,Int64},Adj::SparseMatrixCSC{Int64,Int64}, t::Symbol )
	degree = sum(Adj,dims = 2)[:]
	if t == :Fun
		return collect([EnergyFun((x::Vector{Int64}-> -( J[i,j]*x[1]*x[2]+x[1]*(h[i]/degree[i])+ x[2]*(h[j]/degree[j]))),[i,j] ) for i=1:N, j=1:N if Adj[i,j]!=0 && i<j])
	elseif t == :PairIsing
		return collect([PairIsing(J[i,j],h[i]/degree[i],h[j]/degree[j],[i;j] ) for i=1:N, j=1:N if Adj[i,j]!=0 && i<j])
	end
end

struct DCState
    μtl :: Vector{Vector}			# array of first moments of tilted distributions (one for each factor)
    Σtl :: Vector{Matrix}			# array of second moments of tilted distributions (one for each factor)
    μt  :: Vector					# mean vector w.r.t tilted distribution (size N)
    Σt  :: Matrix					# correlation matrix w.r.t tilted distribution (size NxN)
	µ   :: Vector					# mean vector of gaussian distribution (size N)
    Σ   :: Matrix					# correlation matrix of the gaussian distribution  (size NxN)
    h   :: Vector{Vector}			# array of gaussian field vector (one for each factor)
    S   :: Vector{Matrix}			# array of inverse correlation matrix (one for each factor)
    yc  :: Vector{Vector}			# array of cavity fields vector (one for each factor)
    Sc  :: Vector{Matrix}			# array of cavity inverse correlation matrix (one for each factor)
end

function DCState(N::Int,m::Vector{Int};f::Function=n->2e-2*(rand(n) .- 0.5))
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
