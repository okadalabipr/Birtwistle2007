# Birtwistle_et_al_2007
Birtwistle, M. R. *et al.* Ligand-dependent responses of the ErbB signaling network: Experimental and modeling analyses. *Mol. Syst. Biol.* **3**, (2007). https://doi.org/10.1038/msb4100188

## Requirements
- **[Julia 1.0+](https://julialang.org)**
    - [Sundials](https://github.com/JuliaDiffEq/Sundials.jl)
    - [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
    - [IJulia](https://github.com/JuliaLang/IJulia.jl)

## Run Simulation and View Results
```julia
include("Birtwistle_et_al_2007.jl")
using .Birtwistle_et_al_2007
runSim()
```
![ErbBmodel](https://user-images.githubusercontent.com/31299606/60935404-c150b780-a304-11e9-9c67-8a14e8ba62d7.png)

## Installation

    $ git clone https://github.com/okadalabipr/Birtwistle_et_al_2007.git
