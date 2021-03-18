using JLD, SparseArrays, Test, Random

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
        #βvec = sort(collect(keys(true_moments)))
        βvec = [0.1,0.3]
        @testset "β $β" for β in βvec
            psi = DC.IsingFG(N,β*model["h"],β*sparse(model["J"]), :Fun);
            ϕ , res, it, epsv = DC.density_consistency(psi,N,verbose = false, epsconv = 1e-10,ρ = 0.9)
            @test isapprox(true_moments[β][1], ϕ.µt, atol = 1e-7)
        end
    end
end


# USELESS
#@testset "Reproducibility with rng seeds" begin
#    dir = "test/random_graphs/"
#    ls = readdir(dir)
#    filenames = replace.(ls[findall(occursin.("model",ls))], "model.jld"=>"")
#    model = load(string(dir,filenames[1],"model.jld"));
#    N = length(model["h"])
#    β = 0.2
#    psi = DC.IsingFactors(N,β*model["h"],β*sparse(model["J"]), :Fun);
#
#    Random.seed!(1);
#    ϕstart = DC.DCState(N,[length(psi[i].idx) for i in 1:length(psi)]; f = x -> rand(x) .- 0.5);
#    ϕ1a , res, it, epsv = DC.density_consistency(psi,N,verbose = false, epsconv = 1e-10, seed = 1,state = ϕstart )
#
#    Random.seed!(1);
#   ϕstart = DC.DCState(N,[length(psi[i].idx) for i in 1:length(psi)]; f = x -> rand(x) .- 0.5);
#    ϕ1b , res, it, epsv = DC.density_consistency(psi,N,verbose = false, epsconv = 1e-10, seed = 1,state = ϕstart )
#
#    #ϕ2 , res, it, epsv = DC.density_consistency(psi,N,verbose = false, epsconv = 1e-10, seed = 2)
#    @testset "fields $x" for x in fieldnames(typeof(ϕ1a))
#        @test getfield(ϕ1a,x) == getfield(ϕ1b,x)
#    end
#end


#@testset "Factor Constructors"  begin
#    dir = "test/random_graphs/"
#    ls = readdir(dir)
#    filenames = replace.(ls[findall(occursin.("model",ls))], "model.jld"=>"")
#    @testset "Ising open-chains $i" for i in 1:length(filenames)
#        model = load(string(dir,filenames[i],"model.jld"));
#        N = length(model["h"])
#        βvec = [0.1]
#        @testset "β $β" for β in βvec
#            psiF = DC.IsingFactors(N,β*model["h"],β*sparse(model["J"]), model["Adj"], :Fun);
#            psiPI = DC.IsingFactors(N,β*model["h"],β*sparse(model["J"]), model["Adj"], :PairIsing);
#            ϕF , res, it, epsv = DC.density_consistency(psiF,N,verbose = false, epsconv = 1e-10)
#            ϕI , res, it, epsv = DC.density_consistency(psiPI,N,verbose = false, epsconv = 1e-10)
#            @testset "fields $x" for x in fieldnames(typeof(ϕF))
#                @test isapprox(getfield(ϕF,x),getfield(ϕI,x),atol = 1e-9)
#            end
#        end
#    end
#end
@testset "Equivalence of Ising Constructors" begin
    dir = "test/random_graphs/"
    ls = readdir(dir)
    filenames = replace.(ls[findall(occursin.("model",ls))], "model.jld"=>"")
    model = load(string(dir,filenames[3],"model.jld"));
    N = length(model["h"])
    β = 0.2
    psiF = DC.IsingFG(N,β*model["h"],β*sparse(model["J"]), :Fun);
    psiPI = DC.IsingFG(N,β*model["h"],β*sparse(model["J"]), :PairIsing);
    ϕF , res, it, epsv = DC.density_consistency(psiF,N,verbose = false, epsconv = 1e-10, ρ = 0.9)
    ϕI , res, it, epsv = DC.density_consistency(psiPI,N,verbose = false, epsconv = 1e-10, ρ = 0.9)
    #@testset "fields $x" for x in fieldnames(typeof(ϕF))
    #    @test isapprox(getfield(ϕF,x),getfield(ϕI,x),atol = 1e-9)
    #end
    @test isapprox(ϕF.µt,ϕI.µt,atol = 1e-7)
end
