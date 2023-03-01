import numpy as np
import random

def pso_objective(x):
    F = np.zeros(5)
    # Incorporating all equations and constraints
    # x(1) = SiO2
    # x(2) = O2
    # x(3) = SiO
    # x(4) = H2O
    # x(5) = H2
    # x(6) = FeO
    # x(7) = Fe


    # SiO2 = SiO + .5O2
    F[0] = ((x[0] * (fo2_0**0.5)) / x[1] - 10**log_K_SiO2_picked)**2

    # Si Mass Balance
    F[1] = (1-((4*np.pi*(x[0]*1e5)*R_2/g)*28/44*1000/mag + 
             (x[1]*bse_grams)*28/60/mag) / 
          (sio2_initial*28/60/mag))**2

    F[2] = (1-(x[1]/(x[2]*x_h2o[0]**0.5))/10**log_K_H2O_picked)
    F[3] = (1-((4*np.pi*(x[2]*1e5)*R_2/g)*1000/mag + (4*np.pi*(x[1]*1e5)*R_2/g)*1000*2/18/mag
               + 1e-8*x_h2o[2]*np.exp(-4000/2000)*bse_grams/mag + (bse_grams*0.01*(x[1]/10/241.5)**0.74)*2/18/mag)
             /(H2_i/mag+(H2O_i*2/18)/mag))
    F[4] = 1e3*(1-((4*np.pi*x[1]*1e5*R_2/g)*16/18*1000/mag
                  + (bse_grams*0.01*(x[1]/10/241.5)**0.74)*16/18/mag
                  + (4*np.pi*x[0]*1e5*R_2/g)*1000/mag)
             /(O_bulk_pick/mag))
    sum_squared_error = np.sum(F**2)
    return sum_squared_error

def particle_swarm_optimization(pso_objective, x_min, x_max, num_particles, max_iterations):
    x = np.zeros((num_particles, 3))
    v = np.zeros((num_particles, 3))
    p = np.zeros((num_particles, 3))
    fp = np.zeros(num_particles)
    fg = np.inf
    g = np.zeros(3)
    for i in range(num_particles):
        x[i,:] = x_min + (x_max - x_min) * np.random.rand(3)
        v[i,:] = np.zeros(3)
        p[i,:] = x[i,:]
        fp[i] = pso_objective(x[i,:])
        if fp[i] < fg:
            fg = fp[i]
            g = x[i,:]
    w = 0.9
    c1 = 2.0
    c2 = 2.0
    for iteration in range(max_iterations):
        for i in range(num_particles):
            r1 = np.random.rand(3)
            r2 = np.random.rand(3)
            v[i,:] = w * v[i,:] + c1 * r1 * (p[i,:] - x[i,:]) + c2 * r2 * (g - x[
