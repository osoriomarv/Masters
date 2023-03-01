import numpy as np

def KD_solver(Pressures, Temp):
    a_i_Ni = 0.46  #unitless
    b_i_Ni = 2700 # K
    c_i_Ni = -61 #K/ GPa

    a_i_Si = 1.3 #unitless
    b_i_Si = -13500 #K/GPa

    a_i_O = .6 #Unitless
    b_i_O =-3800 #K
    c_i_O =22 #K/GPa

    log_K_O_D = np.zeros(100)
    log_K_Si_D = np.zeros(100)
    log_K_Ni_D = np.zeros(100)

    for j in range(100):
        log_K_O_D[j] = a_i_O + (b_i_O/Temp[j]) + (c_i_O*Pressures[j])/(Temp[j])
        log_K_Si_D[j] = a_i_Si + (b_i_Si/Temp[j])
        log_K_Ni_D[j] = a_i_Ni + (b_i_Ni/Temp[j]) + (c_i_Ni*Pressures[j])/(Temp[j])

    return log_K_O_D, log_K_Si_D, log_K_Ni_D

log_K_O_D, log_K_Si_D, log_K_Ni_D = KD_solver(Pressures, Temp)
