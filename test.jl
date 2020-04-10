include("src/GEP.jl")
using JLD, SparseArrays

tree = load("test_tree.jld");
psi = DC.IsingConstructor(20,tree["h"],tree["J"],tree["degreevector"],:Fun);
psiI = DC.IsingConstructor(20,tree["h"],tree["J"],tree["degreevector"],:Fun);
ϕ , res, epsv = @timed DC.density_consistency(psi,20)
#ϕI , res, epsv = DC.density_consistency(psiI,20)
#println([tree["µtrace"] ϕ.µt ϕI.µt])
