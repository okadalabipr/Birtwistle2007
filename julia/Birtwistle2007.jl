module Birtwistle2007

using PyPlot

export runSim

include("model/simulation.jl")
using .Sim

include("plotFunc.jl");

function runSim()
    plotFunc_timecourse(Sim)
end

end  # module