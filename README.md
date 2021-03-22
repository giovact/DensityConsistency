# Density Consistency

A julia code to run Density Consistency algorithm.

DC allows to compute marginal marginal distributions of discrete graphical models:

<img src="https://render.githubusercontent.com/render/math?math=p\left(\boldsymbol{x}\right) = \frac{1}{Z}\prod_{a}\psi_{a}\left(\boldsymbol{x}_{a}\right)\prod_{i}\Delta_{i}\left(x_{i}\right)">

where xᵢ ∈ {-1,1}

Each factor ψₐ is defined as ψₐ ∝ exp(-E(xₐ))   a structure with two fields:
* `E::Function`: E(xₐ)
* `idx::Vector{Int64}`: ∂a indeces

Check the script `types.jl` for basic Ising Factor Graph constructors.

### Package requirements
- ExtractMacro
- Random
- SparseArrays

## Usage
The package is not registered. To use it, clone the repository locally:
> git clone https://github.com/giovact/DensityConsistency

### Loading
```julia
include("src/DC.jl")
using .DC
```
A small number of tests can be run by typing `julia run_test.jl` (requires Tests.jl package)
Check `example_Ising.ipynb` jupyter notebook for basic usage on the Ising model

The algorithm runs by calling the function ``density_consistency``, taking the following arguments:

* `Ψ::Vector{<:Factor}`: Array of factor nodes ψₐ
* `N::Int`: number of spins

Some additional named arguments:

* `closure::Symbol = :DC`: closure protocol
* `η::Float64 = 1.0`: interpolation parameter (valid on :DC closure) -> ``η=0`` gives BP fixed points
* `epsconv::Float64 = 1e-6`: precision convergence
* `update::Symbol = :par`: type of update -> parallel (`:par`) or random sequential (`:seq`)
* `ρ::Float64 = 0.9`: damping for parameters' update (`ρ = 0` means no damping)
* `convtype::Symbol = :params`: convergence criterion on `:moments` or  gaussian factor `:params`

### Reference
Alfredo Braunstein, Giovanni Catania and Luca Dall’Asta
*Loop corrections in spin models through density consistency*
2019, Phys. Rev. Lett. [123020604][papero], arXiv:[1810.10602][paperoarxiv]


[papero]: <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.020604>
[paperoarxiv]: <https://arxiv.org/abs/1810.10602>
[example_notebook]: <https://arxiv.org/abs/1810.10602>
