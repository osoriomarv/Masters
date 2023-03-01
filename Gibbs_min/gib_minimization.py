import numpy as np

T = 1000 # K 
R = 8.314e-3 # J/mol/K

P = 10.0 # atm
Po = 1.0 # atm

species = ['CO', 'H2O', 'CO2', 'H2']

# Heat of formation at 298.15 K

Hf298 = [-110.5271, 
        -241.8264, 
        -393.5224, 
            0.0] # kJ/mol
# Shomate parameters for each species
#           A          B           C          D          E            F          G       H
WB = [[25.56759,  6.096130,     4.054656,  -2.671301,  0.131021, -118.0089, 227.3665,   -110.5271],  # CO
      [30.09200,  6.832514,     6.793435,  -2.534480,  0.082139, -250.8810, 223.3967,   -241.8264],  # H2O
      [24.99735,  55.18696,   -33.69137,    7.948387, -0.136638, -403.6075, 228.2431,   -393.5224],  # CO2
      [33.066178, -11.363417,  11.432816,  -2.772874, -0.158558, -9.980797, 172.707974,    0.0]]     # H2


WB = np.array(WB)
print(WB)
t = T/1000
T_H = np.array([t, t**2 / 2.0, t**3 / 3.0, t**4 / 4.0, -1.0 / t, 1.0, 0.0, -1.0])
T_S = np.array([np.log(t), t,  t**2 / 2.0,  t**3 / 3.0, -1.0 / (2.0 * t**2), 0.0, 1.0, 0.0])

H = np.dot(WB, T_H)        # (H - H_298.15) kJ/mol
S = np.dot(WB, T_S/1000.0)        # abs entropy

Gjo = Hf298 + H - T*S

def func(nj):
    nj = np.array(nj)
    Enj = np.sum(nj);
    Gj =  Gjo / (R * T) + np.log(nj / Enj * P / Po)
    return np.dot(nj, Gj)

Aeq = np.array([[ 1,    0,    1,    0],  # C balance
                [ 1,    1,    2,    0],  # O balance
                [ 0,    2,    0,    2]]) # H balance

# equimolar feed of 1 mol H2O and 1 mol CO
beq = np.array([1,  # mol C fed
                2,  # mol O fed
                2]) # mol H fed

def ec1(nj):
    'conservation of atoms constraint'
    return np.dot(Aeq, nj) - beq

from scipy.optimize import fmin_slsqp

n0 = [0.5, 0.5, 0.5, 0.5]  # initial guesses
N = fmin_slsqp(func, n0, f_eqcons=ec1)
print (N)

yj = N / np.sum(N)
Pj = yj * P

for s, y, p in zip(species, yj, Pj):
    if y < 0 or p < 0:
        continue
    print ('{0:10s}: {1:1.2f} {2:1.2f}'.format(s, y, p))

nuj = np.array([-1, -1, 1, 1])  # stoichiometric coefficients of the reaction
K = np.prod(yj**nuj)
print (K)
