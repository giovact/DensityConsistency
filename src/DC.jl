module DC

using ExtractMacro, LinearAlgebra, SparseArrays, Random
export DCState, density_consistency

include("types.jl")
include("Factor.jl")
include("density_consistency.jl")

end # end module
