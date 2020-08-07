module DC

using ExtractMacro, LinearAlgebra, SparseArrays,Random
export DCState

include("types.jl")
include("Factor.jl")
include("dc_main.jl")

end # end module
