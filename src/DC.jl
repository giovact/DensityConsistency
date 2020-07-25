module DC

using ExtractMacro, LinearAlgebra, SparseArrays,Random
export DCState

include("types.jl")
include("Factor.jl")
include("ep.jl")

end # end module
