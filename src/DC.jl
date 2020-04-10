module DC

using ExtractMacro,
    LinearAlgebra,
    SparseArrays,
    Random

include("types.jl")
include("Factor.jl")
include("ep.jl")

end # end module
