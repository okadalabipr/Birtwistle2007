import numpy as np
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import diffeq, param_values, initial_values

def solveode(diffeq, y0, tspan, args):
    sol = ode(diffeq)
    sol.set_integrator(
        'vode', method='bdf', with_jacobian=True,
        atol=1e-9, rtol=1e-9, min_step=1e-8
    )
    sol.set_initial_value(y0, tspan[0])
    sol.set_f_params(args)

    T = [tspan[0]]
    Y = [y0]

    while sol.successful() and sol.t < tspan[-1]:
        sol.integrate(sol.t+1.)
        T.append(sol.t)
        Y.append(sol.y)

    return np.array(T), np.array(Y)


class Simulation(object):

    tspan = range(1801)
    t = np.array(tspan)
    conditions = [
        'EGF00_HRG05',
        'EGF00_HRG10',
        'EGF05_HRG00',
        'EGF05_HRG05',
        'EGF05_HRG10',
        'EGF10_HRG00',
        'EGF10_HRG05',
        'EGF10_HRG10',
    ]
    
    ERK_act = np.empty((len(t), len(conditions)))
    Akt_act = np.empty((len(t), len(conditions)))
    
    x = param_values()
    y0 = initial_values()
    
    for i, condition in enumerate(conditions):
        if condition == 'EGF00_HRG05':
            y0[V.E] = 0.0
            y0[V.H] = 0.5
        elif condition == 'EGF00_HRG10':
            y0[V.E] = 0.0
            y0[V.H] = 10.0
        elif condition == 'EGF05_HRG00':
            y0[V.E] = 0.5
            y0[V.H] = 0.0
        elif condition == 'EGF05_HRG05':
            y0[V.E] = 0.5
            y0[V.H] = 0.5
        elif condition == 'EGF05_HRG10':
            y0[V.E] = 0.5
            y0[V.H] = 10.0
        elif condition == 'EGF10_HRG00':
            y0[V.E] = 10.0
            y0[V.H] = 0.0
        elif condition == 'EGF10_HRG05':
            y0[V.E] = 10.0
            y0[V.H] = 0.5
        elif condition == 'EGF10_HRG10':
            y0[V.E] = 10.0
            y0[V.H] = 10.0
        
        (T, Y) = solveode(diffeq, y0, tspan, tuple(x))
        
        ERK_act[:,i] = Y[:,V.ERKstar] + Y[:,V.pERK_ERKPpase]
        Akt_act[:,i] = Y[:,V.Aktstar]
            