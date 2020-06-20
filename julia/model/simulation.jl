module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")

using .C
using .V
using Sundials
# using ODE

p = param_values()
u0 = initial_values()

const tspan = (0.0,1800.0)
const t = collect(tspan[1]:1.0:tspan[end])

const condition = 8

ERK_act = zeros(length(t),condition)
Akt_act = zeros(length(t),condition)

for i=1:condition
    if i==1
        u0[V.E] = 0.0
        u0[V.H] = 0.5
    elseif i==2
        u0[V.E] = 0.0
        u0[V.H] = 10.0
    elseif i==3
        u0[V.E] = 0.5
        u0[V.H] = 0.0
    elseif i==4
        u0[V.E] = 0.5
        u0[V.H] = 0.5
    elseif i==5
        u0[V.E] = 0.5
        u0[V.H] = 10.0
    elseif i==6
        u0[V.E] = 10.0
        u0[V.H] = 0.0
    elseif i==7
        u0[V.E] = 10.0
        u0[V.H] = 0.5
    elseif i==8
        u0[V.E] = 10.0
        u0[V.H] = 10.0
    end

    prob = ODEProblem(diffeq,u0,tspan,p)
    sol = solve(prob,CVODE_BDF(),saveat=1.0,abstol=1e-9,reltol=1e-9)
    #(sol_t,sol_u) = ode23s(diffeq,u0,t;points=:specified)

    for j=1:length(t)
        ERK_act[j,i] = sol.u[j][V.ERKstar] + sol.u[j][V.pERK_ERKPpase]
        Akt_act[j,i] = sol.u[j][V.Aktstar]
    end

end
end  # module