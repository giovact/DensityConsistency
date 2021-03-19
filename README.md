# Density Consistency

A julia code to run Density Consistency algorithm. 

Check out our paper for additional details!

Alfredo Braunstein, Giovanni Catania and Luca Dall’Asta
*Loop corrections in spin models through density consistency*
2019, Phys. Rev. Lett. [123020604][papero], arXiv:[1810.10602][paperoarxiv]

### Theory
The algorithm allows to compute marginal marginal distributions of:
<img src="https://render.githubusercontent.com/render/math?math=p\left(\boldsymbol{x}\right) = \frac{1}{Z}\prod_{a}\psi_{a}\left(\boldsymbol{x}_{a}\right)\prod_{i}\Delta_{i}\left(x_{i}\right)">
where <img src="https://render.githubusercontent.com/render/math?math=x_{i} \in \{-1,1\}">
### Package requirements
- ExtractMacro
- SparseArrays
- Random

## Usage
The package is not yet registered. To use it, clone the repository locally:

> git clone https://github.com/giovact/DensityConsistency

### Loading
```julia
include("src/DC.jl")
using .DC
```
A small number of tests can be run by typing `julia run_test.jl`
Check the  `example.ipynb` jupyter notebook for basic usage

[papero]: <https://journals-aps-org.ezproxy.biblio.polito.it/prl/abstract/10.1103/PhysRevLett.123.020604>
[paperoarxiv]: <https://arxiv.org/abs/1810.10602>
[example_notebook]: <https://arxiv.org/abs/1810.10602>


