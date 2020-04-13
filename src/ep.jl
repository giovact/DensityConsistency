eye(N) = Matrix(1.0I, N, N)

function update!(old,new,ρ=0.0)
    r=norm(new-old,Inf);
	old .*=ρ ;
	old .+= (1-ρ)*new;
    return r
end

"""
Density consistency or generic consistency propagation algorithm
"""

function density_consistency(Ψ::Vector{<:Factor},N::Int;
			algvars::DCparams = DCparams(:DC),
			Hg::Vector{Float64} = zeros(N),
			Ag::Matrix{Float64} = zeros(N,N),
			convtype::Symbol = :moments,
			callback::Function  = (x...)->nothing,
            state::DCState = DCState( N, [length(Ψ[a].idx) for a=1:length(Ψ)]) )

    @extract state : μtl Σtl μt Σt µ Σ h S yc Sc
	@extract algvars : maxiter η γ0 epsconv Λ λ closure update rndamp epsclamp ρ verbose

	A,y = zeros(N,N), zeros(N)
    Id = Matrix(1.0I, N, N); M = length(Ψ);
    ε,εₘ = 0.0 , 0.0
    iter = 1
    γ = 0.0
    while iter<=maxiter
    	γ += γ0
        A[:] .= 0.0; y .= μt * γ #check µt?
        A .+= Ag; y .+= Hg
        μt[:] .= 0.0; Σt[:] .= 0.0; #may be useless
        ε,εₘ = 0.0 ,  0.0
        for a in 1:M
            ∂a = Ψ[a].idx #might be put before the while
            A[∂a,∂a] .+= S[a]; y[∂a] .+= h[a]
        end
        perm = collect(1:M) #might be put at the beginning
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
            ynew, Snew, εₘ = setclosure!(ψₐ,yc[a],Sc[a],µtl[a],Σtl[a],closure,ρ,λ,epsclamp,εₘ)
            if update == :seq
                ε = max(ε, update!(S[a], Snew, η), update!(h[a], ynew, η))
				A[∂a,∂a] .+= S[a]; y[∂a] .+= h[a]
            elseif update == :par
                ηr = rndamp*(1-ρ)*(2*rand()-1) + η    # ρr ~ U(2ρ-1,1) (with <ρ̃>=ρ) if rndamp true
                ε = max(ε, update!(S[a], Snew,ηr), update!(h[a], ynew,ηr))
            end
            μt[∂a] .= µtl[a] ; Σt[∂a,∂a] .= Σtl[a]
        end
        callback(state,iter,[ε, εₘ])
        iter+=1
        if check_convergence(ε,εₘ,epsconv,convtype)
			verbose && println("Closure $closure converged in $iter iterations w.r.t $convtype")
			return state, :converged, iter, [ε, εₘ]
		end
    end
	verbose && println("Closure $closure NOT converged w.r.t $convtype")
    return state, :unconverged, iter, [ε, εₘ]
end

function check_convergence(err_params::Float64,err_moments::Float64,epsconv::Float64,convtype::Symbol)
	convtype == :params && return err_params < epsconv
	convtype == :moments && return err_moments < epsconv
end
