module DC

using ExtractMacro, LinearAlgebra, SparseArrays,Random
export DCState, density_consistency

include("types.jl")
include("Factor.jl")
include("dc_main.jl")

end # end module
