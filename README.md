# Birtwistle_et_al_2007
Birtwistle, M. R. *et al.* Ligand-dependent responses of the ErbB signaling network: Experimental and modeling analyses. *Mol. Syst. Biol.* **3**, (2007). https://doi.org/10.1038/msb4100188

## Requirements
- **[Julia 1.0+](https://julialang.org)**
  - [ODE](https://github.com/JuliaDiffEq/ODE.jl)
  - [PyPlot](https://github.com/JuliaPy/PyPlot.jl)
- **[Juno](http://junolab.org)**

## Run Simulation and View Results
```julia
include("simulation.jl");
include("plot_func.jl");
savefig("./ErbBmodel.png",bbox_inches="tight");
```
## Installation

    $ git clone https://github.com/okadalabipr/Birtwistle_et_al_2007.git
