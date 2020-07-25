using JLD, SparseArrays, Test

include("src/DC.jl")
using .DC

@testset "exactness on trees" begin
    tree = load("test/ising_chain.jld");
    N = length(tree["h"])
    psi = DC.IsingConstructor(20,tree["h"],sparse(tree["J"]),:Fun);
    ϕ , res, it, epsv = DC.density_consistency(psi,N)
    @test isapprox(tree["mean"], ϕ.µt, atol = 1e-6)
end
#ϕI , res, epsv = DC.density_consistency(psiI,20)
#println([tree["µtrace"] ϕ.µt ϕI.µt])
