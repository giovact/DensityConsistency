abstract type Factor end

struct FactorFun <: Factor
    E::Function
    idx::Vector{Int64}
end

struct FactorIsing <: Factor
    J::Float64
    hi::Float64
    hj::Float64
    idx::Vector{Int64}
end

function IsingConstructor(N::Int64,h::Vector{Float64},J::Matrix{Float64},d::Vector{Int64},t::Symbol)
	if t == :Fun
		return collect([FactorFun((x::Vector{Int64}-> -( J[i,j]*x[1]*x[2]+x[1]*(h[i]/d[i])+ x[2]*(h[j]/d[j]))),[i,j] ) for i=1:N, j=1:N if J[i,j]!=0 && i<j])
	elseif t== :Is
		return collect([FactorIsing(J[i,j],h[i]/d[i],h[j]/d[j],[i;j] ) for i=1:N, j=1:N if J[i,j]!=0 && i<j])
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

function DCState(N::Int,m::Vector{Int};f::Function=n->(rand(n) .- 0.5))
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


struct DCparams
	maxiter	:: Int				# maximum number of iterations
	ρ :: Float64				# damping for parameters' update (ρ = 0 means no damping)
	γ0:: Float64				# reinforcement
	epsconv :: Float64			# precision convergence
	Λ :: Float64				# add a diagonal matrix Λ*I in the inversion of the full correlation matrix
	λ :: Float64				# add a diagonal matrix λ*I in the inversion of the cavity and tilted correlation matrices
	closure :: Symbol			# closure protocol (default = :DC)
	update :: Symbol			# type of update (parallel or sequential)
	rndamp :: Bool				# eventually, apply a random damping (might be useful in parallel update)
	epsclamp :: Float64			# clamp
	η :: Float64				# interpolation parameter (DC closure)
	verbose :: Bool				# print at (un)convergence
end


function DCparams(closure::Symbol;maxiter::Int64 = 2000,
					ρ::Float64 = 0.9,
					γ0::Float64 = 0.0,
					epsconv::Float64 = 1e-6,
					Λ::Float64 = 0.0,
					λ::Float64 = 0.0,
					update::Symbol = :par,
					rndamp::Bool = false,
					epsclamp::Float64 = 0.0,
					η::Float64 = 1.0,
					verbose::Bool = true)

	return DCparams(maxiter,ρ,γ0,epsconv,Λ,λ,closure,update,rndamp,epsclamp,η,verbose)
end
