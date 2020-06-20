module Birtwistle2007

using PyPlot;

export runSim

include("model/name2idx/parameters.jl")
include("model/name2idx/species.jl")
include("model/set_model.jl")
include("model/simulation.jl")

using .Sim

include("plotFunc.jl");

function runSim()
    plotFunc_timecourse(Sim)
end

end  # module