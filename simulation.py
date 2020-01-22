import numpy as np
from scipy.integrate import ode

from model.name2idx import parameters as C
from model.name2idx import variables as V
from model.param_const import f_params
from model.initial_condition import initial_values
from model.differential_equation import diffeq


def solveode(diffeq,y0,tspan,args):
    sol = ode(diffeq)
    sol.set_integrator('vode',method='bdf',min_step=1e-8,with_jacobian=True)
    sol.set_initial_value(y0,tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < tspan[-1]:
        sol.integrate(sol.t+1.)
        T.append(sol.t)
        Y.append(sol.y)

    return np.array(T),np.array(Y)


class Simulation(object):

    tspan = range(1801)
    condition = 8
    
    t = np.array(tspan)
        
    ERK_act = np.empty((len(t),condition))
    Akt_act = np.empty((len(t),condition))
    
    x = f_params()
    y0 = initial_values()
    
    for i in range(condition):
        if i==0:
            y0[V.E] = 0.0
            y0[V.H] = 0.5
        elif i==1:
            y0[V.E] = 0.0
            y0[V.H] = 10.0
        elif i==2:
            y0[V.E] = 0.5
            y0[V.H] = 0.0
        elif i==3:
            y0[V.E] = 0.5
            y0[V.H] = 0.5
        elif i==4:
            y0[V.E] = 0.5
            y0[V.H] = 10.0
        elif i==5:
            y0[V.E] = 10.0
            y0[V.H] = 0.0
        elif i==6:
            y0[V.E] = 10.0
            y0[V.H] = 0.5
        elif i==7:
            y0[V.E] = 10.0
            y0[V.H] = 10.0
        
        (T,Y) = solveode(diffeq,y0,tspan,tuple(x))
        
        ERK_act[:,i] = Y[:,V.ERKstar] + Y[:,V.pERK_ERKPpase]
        Akt_act[:,i] = Y[:,V.Aktstar]
            