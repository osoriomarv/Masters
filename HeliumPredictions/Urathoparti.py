import numpy as np
import matplotlib.pyplot as plt

U_amu = 238.02891
U_concentration_mol_g = 0.02 / 1e6 / U_amu
U_238_mol_g = U_concentration_mol_g * 0.99275
U_235_mol_g = U_concentration_mol_g * 0.00720

Th_amu = 232.03806
Th_232_mol_g = 0.08 / 1e6 / Th_amu

U_238_lambda = -1.54e-10
U_235_lambda = -9.72e-10
Th_232_lambda = -4.940e-11

U_238_initial_mol_g = U_238_mol_g / np.exp(-U_238_lambda * 4.5e9)
U_235_initial_mol_g = U_235_mol_g / np.exp(-U_235_lambda * 4.5e9)
Th_232_initial_mol_g = Th_232_mol_g / np.exp(-Th_232_lambda * 4.5e9)

atm_3_4_he = 1.34e-6

delIW = -6

T_adiabat = 1833.3777 + (((4500 * (1 / 4)) * 1))

X_sio2 = 0.456
X_s = 0.0355

a_coeff = 0.1
b_coeff = -11000
c_coeff = -4700
d_coeff = 13000
e_coeff = -0.49

BSE_initial = 1.68e-10
BSE_initial_Th = 6.90e-10

logD = a_coeff + (b_coeff / T_adiabat) + c_coeff * ((1 - X_s) ** 2) / T_adiabat + d_coeff * ((1 - X_sio2) ** 2) / T_adiabat + e_coeff * delIW

Dm = 10 ** logD

He_ratio = 1.66e-4

U_mol_g_metal = (BSE_initial * Dm)
U_mol_metal = U_mol_g_metal * 1.9831e+27
U235_mol_metal = U_mol_metal * 0.00720
U238_mol_metal = U_mol_metal * 0.99275

Th_mol_g_metal = BSE_initial_Th * Dm
Th_mol_metal = Th_mol_g_metal * 1.9831e+27

he_4 = 1.69084641083620e+17 / 4
he_3 = He_ratio * he_4

time_vector = np.arange(0, 4.5e9, 1e6)

he3_ratio_temporal = []

for t in time_vector:
    U238_left = U238_mol_metal * np.exp(-1.54e-10 * t)
    U235_left = U235_mol_metal * np.exp(-9.72e-10 * t)
    Th232_left = Th_mol_metal * np.exp(-4.940e-11 * t)

    alpha_Th232_mol = (Th_mol_metal - Th232_left) * (6)
    alpha_U235_mol = (U235_mol_metal - U235_left) * (7)
    alpha_U238_mol = (U238_mol_metal - U238_left) * (8)

    He_4_temporal = he_4 + alpha_Th232_mol + alpha_U235_mol + alpha_U238_mol
    he3_ratio_temporal.append(he_3 / He_4_temporal)

he3_ratio_temporal = np.array(he3_ratio_temporal)

plt.plot(time_vector, (he3_ratio_temporal / atm_3_4_he))
plt.xlabel('Time (GA)')
plt.ylabel('3He/4He')
plt.show()

x = 1e-6 * (1.9831e+24)
