include("src/GEP.jl")
using JLD, SparseArrays

tree = load("test_tree.jld");
psi = DC.IsingConstructor(20,tree["h"],tree["J"],tree["degreevector"],:Fun);
ϕ , res, epsv = DC.density_consistency(psi,20)
