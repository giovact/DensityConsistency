module DensityConsistency

using ExtractMacro, LinearAlgebra, SparseArrays, Random
export DCState, density_consistency

include("types.jl")
include("closures.jl")
include("dc_main.jl")

end # end module
