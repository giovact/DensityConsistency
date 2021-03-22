eye(N) = Matrix(1.0I, N, N)

function update!(old,new,ρ=0.0)
    r=norm(new-old,Inf);
	old .*=ρ ;
	old .+= (1-ρ)*new;
    return r
end


"""
Density consistency algorithm

Compute approximate marginals of
``p( \\bf{x} )=\\frac1Z \\Phi(\\bf{x} ) \\prod_a \\psi_{a}(x_{\\partial a})``

where `\\Phi(\\bf{x} )` is an optional multivariate Gaussian distribution

Input:

* `Ψ::Vector{<:Factor}`: Array of factor nodes ``\\psi_{a}``
* `N::Int`: number of spins

Optional named arguments:

* `closure::Symbol = :DC`: closure protocol
* `maxiter::Int64 = 2000`: max number of iterations
* `η::Float64 = 1.0`: interpolation parameter (valid on :DC closure) -> ``η=0`` gives BP fixed points
* `γ0::Float64 = 0.0`: reinforcement -> at each iteration t adds a field vector ``h^{t}_{i} = \\gamma_0 t \\mu^{t}_{i}``
* `epsconv::Float64 = 1e-6`: precision convergence
* `Λ::Float64 = 1e-15`: adds a diagonal matrix ``\\Lambda \\mathbb{I}`` in the inversion of the full correlation matrix
* `λ::Float64 = 1e-15`: adds a diagonal matrix ``\\lambda \\mathbb{I}`` in the inversion of the cavity and tilted correlation matrices
* `update::Symbol = :par`: type of update -> parallel (:par) or random sequential (:seq)
* `rndamp::Bool = false`: a random damping (might be useful in parallel update)
* `epsclamp::Float64 = 1e-15:` clamp tilted moments in [-1+epsclamp, 1-epsclamp]
* `ρ::Float64 = 0.9`: damping for parameters' update (ρ = 0 means no damping)
* `verbose::Bool = true` print at (un)convergence
* `Hg::Vector{Float64} = zeros(N)`: field vector of `\\Phi(\\bf{x} )`
* `Ag::Matrix{Float64} = zeros(N,N)`: inverse covariance matrix of `\\Phi(\\bf{x} )`
* `convtype::Symbol = :params`: convergence criterion on :moments or  gaussian factor :params
* `seed :: Int64 = -1`: seed for RNG
* `callback::Function  = (x...)->nothing`: callback function to print progress report
* `h0scale::Float64 = 1e-3`: field scale for gaussian parameter initialization
* `state::DCState = DCState( N, [length(Ψ[a].idx) for a=1:length(Ψ)]; h0 = h0scale)`: initial state of DC run

Output:

* `state::DCState`: final DCstate -> see [`DCState`](@ref)
* `convergence_flag::Symbol`: :converged or :unconverged with criterion `convtype`
* `iter::Int`: number of iterations needed to converge
* `errors::Vector{Float64}`: errors on params and moments at convergence `[error_params , error_moments]`

"""
function density_consistency(Ψ::Vector{<:Factor},
			N::Int;
			closure::Symbol = :DC,
			maxiter::Int64 = 2000,
			η::Float64 = 1.0,
			γ0::Float64 = 0.0,
			epsconv::Float64 = 1e-6,
			Λ::Float64 = 1e-15,
			λ::Float64 = 1e-15,
			update::Symbol = :par,
			rndamp::Bool = false,
			epsclamp::Float64 = 1e-15,
			ρ::Float64 = 0.9,
			verbose::Bool = true,
			Hg::Vector{Float64} = zeros(N),
			Ag::Matrix{Float64} = zeros(N,N),
			convtype::Symbol = :params,
			seed :: Int64 = -1,
			callback::Function  = (x...)->nothing,
			h0scale::Float64 = 1e-3,
			state::DCState = DCState( N, [length(Ψ[a].idx) for a=1:length(Ψ)]; h0 = h0scale),
			)

    @extract state : μtl Σtl μt Σt µ Σ h S yc Sc

	seed != -1 && Random.seed!(seed) #set the seed

	A,y = zeros(N,N), zeros(N)
    Id = Matrix(1.0I, N, N); M = length(Ψ);
    ε,εₘ = 0.0 , 0.0
    iter = 1
    γ = 0.0
    while iter<=maxiter
    	γ += γ0
        A[:] .= 0.0; y .= μt * γ
        A .+= Ag; y .+= Hg
        μt[:] .= 0.0; Σt[:] .= 0.0;
        ε,εₘ = 0.0 ,  0.0
        for a in 1:M
            ∂a = Ψ[a].idx
            A[∂a,∂a] .+= S[a]; y[∂a] .+= h[a]
        end
        perm = collect(1:M)
        randperm!(perm)
      	Σ .= inv(A + Λ * I); μ .= Σ * y; # also defined in rndseq update
        for a in perm
            ψₐ=Ψ[a]; ∂a = ψₐ.idx
	    	if update == :seq
	    		A[∂a, ∂a] .-= S[a]; y[∂a] .-= h[a]
            	Σa = ((A + Λ * I) \ Id[:, ∂a])[∂a,:]; Sc[a] .= (Σa + λ * I) \ eye(length(∂a))
            	μa = ((A + Λ * I) \ y)[∂a]; yc[a] .= Sc[a] * μa
	    	elseif update == :par
                Sc[a] .= (Σ[∂a, ∂a] + λ * I)\eye(length(∂a)) - S[a]
                yc[a] .= (Σ[∂a, ∂a] + λ * I)\μ[∂a] - h[a]
	    	end
            ynew, Snew, εₘ = setclosure!(ψₐ,yc[a],Sc[a],µtl[a],Σtl[a],closure,η,λ,epsclamp,εₘ)
            if update == :seq
                ε = max(ε, update!(S[a], Snew, ρ), update!(h[a], ynew, ρ))
				A[∂a,∂a] .+= S[a]; y[∂a] .+= h[a]
            elseif update == :par
                ρr = rndamp*(1-ρ)*(2*rand()-1) + ρ    # ρr ~ U(2ρ-1,1) (with <ρ̃>=ρ) if rndamp true
                ε = max(ε, update!(S[a], Snew,ρr), update!(h[a], ynew,ρr))
            end
            μt[∂a] .= µtl[a] ; Σt[∂a,∂a] .= Σtl[a]
        end
        callback(state,iter,[ε, εₘ])
        if check_convergence(ε,εₘ,epsconv,convtype)
			verbose && println("Closure $closure (η = $η) converged in $iter iterations w.r.t $convtype - err params = ", ε, " - err moments = ", εₘ)
			return state, :converged, iter, [ε, εₘ]
		end
		iter+=1
    end
	verbose && println("Closure $closure (η = $η) NOT converged w.r.t $convtype - err params = ", ε, " - err moments = ", εₘ)
    return state, :unconverged, iter, [ε, εₘ]
end

function check_convergence(err_params::Float64,err_moments::Float64,epsconv::Float64,convtype::Symbol)
	convtype == :params && return err_params < epsconv
	convtype == :moments && return err_moments < epsconv
end
