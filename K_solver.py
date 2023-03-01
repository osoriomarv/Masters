# Solve for the exchange coefficient of reactions using gibbs free energy minimization 

# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import math 

# Define the temperature range
T = np.linspace(900,1000,10) #K
P=1  #pascal
Po=1 #pascal
R = 8.314 # J/mol/K
T_ref = 298.15 # K
# Define constants constants

#Pressure
        #FeO    Fe     SiO2
Molar_v=[1.2040, 0.7092, 2.2608] #J/bar

# Define the stoichiometric matrix
VdP_FeO = []
VdP_SiO2 = []
VdP_H2O = []
VdP_H2O_SiO2 = []

# Define the pressure based on the temperature
for i in range(len(T)):
    #FeO
    VdP_FeO.append(Molar_v[0]*P - (Molar_v[1]*P))

    #SiO2
    VdP_SiO2.append((.5*(R*T[i]*math.log(P/Po)) + (R*T[i]*math.log(P/Po))) - (Molar_v[2]*P))

    #H2O
    VdP_H2O.append((R*T[i]*math.log(P/Po)) - (.5*(R*T[i]*math.log(P/Po))  + R*T[i]*math.log(P/Po)))

    #H2O + SiO
    VdP_H2O_SiO2.append(R*T[i]*math.log(P/Po) + (Molar_v[2]*P))

# Load the data
#       A       B       C       D       E       F       G       H       
WB = [[35.69893, 1.731252, -0.509348, 0.059404, -1.248055, -114.6019, 249.1911, -100.416],      #SiO
        [20.9111, 10.7207, -2.0205,     0.1464, 9.2457,     5.3377, 237.6185,   0],             #O2
        [85.772, -0.000016, 0.00004, -3.81E-07, -0.000017, -952.87, 113.344,   -902.661],       #SiO2
        [41.96426, 8.622053, -1.49978, 0.098119, -11.15764, -272.1797, 219.7809, -241.8264],    #H2O 
        [43.41356, -4.293079, 1.272428, -0.096876, -20.533862, -38.515158, 162.081354, 0],      #H2
        [46.024, -1.88E-08, 6.09E-09, -6.64E-10, -8.25E-09, -10.80543, 72.54094, 12.39502],     #Fe
        [68.1992, -4.50E-10, 1.20E-10, -1.06E-11, -3.09E-10, -281.4326, 137.8377, -249.5321]]   #FeO
WB = np.array(WB)

# Entries are in the following order:
S_298= np.array([211.58, #SiO
        205.15, #O2
        47.92,  #SiO2
        188.84, #H2O
        130.68, #H2
        34.76,  #Fe
        75.4]) # FeO 
        #J/mol/K

# Enthalpy of formation at 298.15 K
Hf_298= np.array([-100.416, #SiO
            0, #O2
        -902.661, #SiO2
        -241.8264, #H2O
            0, #H2
        12.39502, #Fe
        -249.5321])*1000# FeO
        #J/mol

# Enthalpy of formation at 298.15 K
# FeO = Fe + .5O2
Hrxn_29815_FeO = (Hf_298[6] - Hf_298[5] - .5*Hf_298[1])
Srxn_29815_FeO = (S_298[6] - S_298[5] - .5*S_298[1]) #J/mol/K
Grxn_29815_FeO = Hrxn_29815_FeO - 298.15*(Srxn_29815_FeO)#J/mol
K_29815_FeO = (-Grxn_29815_FeO/(2.2303*R*T_ref)) #J/mol/K

# H2 + .5O2 = H2O
Hrxn_29815_H2O = (Hf_298[3] - (.5*Hf_298[1] + Hf_298[4])) #J/mol
Srxn_29815_H2O = S_298[3] - (.5*S_298[1]+ S_298[4]) #J/Mol
Grxn_29815_H2O = Hrxn_29815_H2O - 298.15*(Srxn_29815_H2O) #J/mol
K_29815_H2O = (-Grxn_29815_H2O/(2.2303*R*T_ref))

#SiO2 = SiO + .5O2
Hrxn_29815_SiO2 = (.5*Hf_298[1] + Hf_298[1]) - Hf_298[2] #J/mol
Srxn_29815_SiO2 = (.5*S_298[1]+ S_298[1]) - S_298[2] #J/Mol
Grxn_29815_SiO2 = Hrxn_29815_SiO2 - 298.15*(Srxn_29815_SiO2) #J/mol
K_29815_SiO2 = (-Grxn_29815_SiO2/(2.2303*R*T_ref))

# Using vectorized shomate equations to solve for enthalpy and entropy
T_H = np.array([T/1000,  (T/1000)**2 / 2.0, (T/1000)**3 / 3.0, (T/1000)**4 / 4.0, -1.0 / (T/1000), 1.0, 0.0, -1.0])*1000
T_S = np.array([np.log(T/1000), T/1000,  (T/1000)**2 / 2.0,  (T/1000)**3 / 3.0, -1.0 / (2.0 * (T/1000)**2), 0.0, 1.0, 0.0])
H = np.dot(WB, T_H)        # (H - H_298.15) J/mol
S = np.dot(WB, T_S)        # absolute entropy J/mol/K

# create stacked matrix of outputs for each species
S_Species = np.vstack((S))
H_Species = np.vstack((H))
print(H[0])


# Solving for enthalpy and entropy of FeO = Fe + .5O2 at different temperatures
Hrxn_FeO = (Hrxn_29815_FeO + (H_Species[6] - H_Species[5] - .5*H_Species[1])) # j/mol
print(Hrxn_FeO)
# delHrxnlogFeO.append(Hrxn_FeO) # Logs the Hrxn FeO calculations
# Srxn_FeO = S_Species[i][7] - S_Species[i][6] - .5*S_Species[i][2] # J/mol
# delSrxnlogFeO.append(Srxn_FeO) # Logs the entropy 
# Grxn_FeO_Jmol.append((Hrxn_FeO - T[i]*(Srxn_FeO) + VdP_FeO[i])) # J/mol
# log_K_FeO.append((-Grxn_FeO_Jmol[i]/(2.303*R*T[i])))


# def enthalpy(T):
#     T_H = np.array([T/1000,  (T/1000)**2 / 2.0, (T/1000)**3 / 3.0, (T/1000)**4 / 4.0, -1.0 / (T/1000), 1.0, 0.0, -1.0])*1000
#     H = np.dot(WB, T_H)        # (H - H_298.15) J/mol
#     return H

# def entropy(T):
#     T_S = np.array([np.log(T/1000), T/1000,  (T/1000)**2 / 2.0,  (T/1000)**3 / 3.0, -1.0 / (2.0 * (T/1000)**2), 0.0, 1.0, 0.0])
#     S = np.dot(WB, T_S)        # absolute entropy J/mol/K
#     return S

# def gibbs(T):
#     G = enthalpy(T) - T*entropy(T)
#     return G




#Temperature and Pressure corrections
# for i in range(len(S)):
#     Hrxn_FeO = (Hrxn_29815_FeO + (H[i,7] -H[i,6] - .5*H[i,2]))
#     Srxn_FeO = S_298[i,7] - S_298[i,6] - .5*S_298[i,2] #J/mol
#     Grxn_FeO_Jmol = (Hrxn_FeO - T*(Srxn_FeO) + VdP_FeO)
#     log_K_FeO = (-Grxn_FeO_Jmol/(2.303*R*T))
#     print(log_K_FeO)


# Srxn_FeO = S_Species(i,7) - S_Species(i,6) - .5*S_Species(i,2); %J/mol
# delSrxnlogFeO(i)=Srxn_FeO; %Logs the entropy 
# Grxn_FeO_Jmol(i)= (Hrxn_FeO - T(i).*(Srxn_FeO) + VdP_FeO(i));% J/mol
# log_K_FeO(i) = (-Grxn_FeO_Jmol(i)./(2.303.*R.*T(i)));

# Gjo = Hf_298 + H - T*S      # Gibbs energy of each component at 1000 K

# print(Gjo)

# def func(nj):
#     nj = np.array(nj)
#     Enj = np.sum(nj);
#     Gj =  Gjo / (R * T) + np.log(nj / Enj * P / Po)
#     return np.dot(nj, Gj)
