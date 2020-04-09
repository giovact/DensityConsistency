eye(N) = Matrix(1.0I, N, N)

function update!(old,new,ρ=0.0)
    r=norm(new-old,Inf); old .*=ρ ; old .+= (1-ρ)*new;
    return r
end

"""
Density consistency or generic consistency propagation algorithm
"""

function density_consistency(Ψ::Vector{<:Factor},N::Int;
			algvars::DCparams = DCparams(:DC),
			Hg::Vector{Float64} = zeros(N),
			Ag::Matrix{Float64} = zeros(N,N),
			callback::Function  = (x...)->nothing,
            state::DCState = DCState( N, [length(Ψ[a].idx) for a=1:length(Ψ)]) )

    @extract state : μtl Σtl μt Σt µ Σ h S yc Sc
	@extract algvars : maxiter ρ γ0 epsconv Λ λ closure update rndamp epsclamp η verbose

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
        randperm!(perm) #requires Random
      	if update == :par Σ .= inv(A + Λ * I); μ .= Σ * y; end #might be useful also for rndseq
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
            avt, covt = moments(ψₐ,yc[a],Sc[a])
            ynew,Snew = setclosure(ψₐ,yc[a],Sc[a],closure,η,λ,epsclamp)
            if update == :seq
                ε = max(ε, update!(S[a], Snew, ρ), update!(h[a], ynew, ρ))
				A[∂a,∂a] .+= S[a]; y[∂a] .+= h[a]
            elseif update == :par
                ρr= rndamp*(1-ρ)*(2*rand()-1) + ρ    # ρr ~ U(2ρ-1,1) (with <ρ̃>=ρ) if rndamp true
                ε = max(ε, update!(S[a], Snew,ρr), update!(h[a], ynew,ρr))
            end
            εₘ = max(εₘ , update!(μtl[a],avt,0.0), update!(Σtl[a],covt,0.0))
            μt[∂a] .= avt ; Σt[∂a,∂a] .= covt
            #μt[∂a] .+= avt./d[∂a]; Σt[∂a,∂a] .= covt
            #Σt[∂a,∂a] .+= Diagonal(covt).*(1 ./d[∂a]-1)
        end
        #if callback(state,iter,[ε, εₘ]) != nothing
        #    break
        #end
		callback(state,iter,[ε, εₘ])
        iter+=1
        if ε < epsconv
			verbose && println("Closure $closure converged in $iter iterations")
			return state, :converged, iter, [ε, εₘ]
		end
    end
	verbose && println("Closure $closure NOT converged ")
    return state, :unconverged, iter, [ε, εₘ]
end
