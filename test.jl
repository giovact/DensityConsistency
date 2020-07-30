using JLD, SparseArrays, Test

include("src/DC.jl")
using .DC

@testset "Exactness on trees" begin
    dir = "test/Ising_trees/"
    ls = readdir(dir)
    filenames = replace.(ls[findall(occursin.("model",ls))], "model.jld"=>"")
    @testset "Ising open-chains $i" for i in 1:length(filenames)
        model = load(string(dir,filenames[i],"model.jld"));
        true_moments = load(string(dir,filenames[i],"moments.jld"),"moments");
        N = length(model["h"])
        βvec = sort(collect(keys(true_moments)))

        @testset "β $β" for β in βvec
            psi = DC.IsingConstructor(N,β*model["h"],β*sparse(model["J"]),:Fun);
            ϕ , res, it, epsv = DC.density_consistency(psi,N,verbose = false, epsconv = 1e-10, convtype = :moments)
            @test isapprox(true_moments[β][1], ϕ.µt, atol = 1e-7)
        end
    end
end
