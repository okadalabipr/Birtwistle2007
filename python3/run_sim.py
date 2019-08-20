from simulation import Simulation
import plot_func

def run_simulation():
    sim = Simulation()
    sim.numerical_integration()

    plot_func.timecourse(sim)


if __name__ == "__main__":
    run_simulation()
    
'''
%matplotlib inline
from run_sim import run_simulation
run_simulation()
'''