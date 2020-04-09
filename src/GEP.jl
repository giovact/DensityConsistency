module DC

using ExtractMacro,
    LinearAlgebra,
    SparseArrays,
    Random


include("types.jl")
include("Factor.jl")
include("ep.jl")

#export FactorFun, FactorIsing

end # end module
