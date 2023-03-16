import numpy as np
import matplotlib.pyplot as plt
import math

Fe_wt = 40.0
Ni_wt = 20.0
Si_wt = 30.0
O_wt = 10.0

# Define atomic masses
atomic_masses = {'Fe': 55.845, 'Ni': 58.6934, 'Si': 28.0855, 'O': 15.999}

total_weight = (Fe_wt + Ni_wt + Si_wt + O_wt) / 100

molar_fracs = {element: (wt / atomic_masses[element]) / total_weight
                for element, wt in zip(atomic_masses, [Fe_wt, Ni_wt, Si_wt, O_wt])}

avg_metal_atomic_mass = sum(molar_fracs[element] * atomic_masses[element] for element in molar_fracs)

SiO2_wt =48.88
TiO2_wt =1.73
Al2O3_wt = 16.77
FeO_wt = 9.71
Fe2O3_wt = 0.00
MgO_wt = 6.65
CaO_wt = 9.86
Na2O_wt = 3.62
K2O_wt = 1.93
# Define atomic masses
atomic_masses = {'Si': 28.085, 'Ti': 47.867, 'Al': 26.982, 'Fe': 55.845, 'O': 15.999, 'Mg': 24.305, 'Ca': 40.078, 'Na': 22.990, 'K': 39.098}

oxide_formulas = {'SiO2': ['Si', 'O', 'O'],
                    'TiO2': ['Ti', 'O', 'O'],
                    'Al2O3': ['Al', 'Al', 'O', 'O', 'O'],
                    'FeO': ['Fe', 'O'],
                    'Fe2O3': ['Fe', 'Fe', 'O', 'O', 'O'],
                    'MgO': ['Mg', 'O'],
                    'CaO': ['Ca', 'O'],
                    'Na2O': ['Na', 'Na', 'O'],
                    'K2O': ['K', 'K', 'O']}

oxide_wts = [SiO2_wt, TiO2_wt, Al2O3_wt, FeO_wt, Fe2O3_wt, MgO_wt, CaO_wt, Na2O_wt, K2O_wt]
total_weight = sum(oxide_wts)

molar_fracs = {oxide: (wt / sum(atomic_masses[element] for element in oxide_formulas[oxide])) / total_weight
                for oxide, wt in zip(oxide_formulas, oxide_wts)}

avg_bse_atomic_mass = sum(molar_fracs[oxide] * sum(atomic_masses[element] for element in oxide_formulas[oxide])
                        for oxide in molar_fracs)


# Calculate the initial concentration of Uranium and thorium
U_amu = 238.02891
Th_amu = 232.03806
# Input parameters
t = 4.64 #Billion years ago

U_concetration_mol_g = 0.02 / 1e6 / U_amu
Th_concentration_mol_g = 0.08 / 1e6 / Th_amu
mass_bse = 3.9663e+27
mass_core = 1.9832e+27

uranium_decay_constant = 4.916e-18  # 238U decay constant (1/yr)
thorium_decay_constant = 1.55e-18  # 232Th decay constant (1/yr)
time_elapsed = 4.5e9  # Time elapsed in years (e.g., age of Earth)

initial_uranium_concentration = U_concetration_mol_g / math.exp(-uranium_decay_constant * time_elapsed)
initial_thorium_concentration = Th_concentration_mol_g / math.exp(-thorium_decay_constant * time_elapsed)

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


logD = a_coeff + (b_coeff / T_adiabat) + c_coeff * ((1 - X_s) ** 2) / T_adiabat + d_coeff * ((1 - X_sio2) ** 2) / T_adiabat + e_coeff * delIW

Dm = 10 ** logD

He_ratio = 1.66e-4

U_mol_g_metal = (initial_uranium_concentration * Dm)
U_mol_metal = U_mol_g_metal * 1.9831e+27
U235_mol_metal = U_mol_metal * 0.00720
U238_mol_metal = U_mol_metal * 0.99275

Th_mol_g_metal = initial_thorium_concentration * Dm
Th_mol_metal = Th_mol_g_metal * 1.9831e+27

he_4 = 1.69084641083620e+17 / 4
he_3 = He_ratio * he_4

time_vector = np.arange(0, 4.5e9, 1e6)

he3_ratio_temporal = []

U_238_mol_g= initial_uranium_concentration*.99275
U_235_mol_g=initial_thorium_concentration*.00720

lambda_U235 = 9.85e-10 # Decay constant of uranium-235 in years^-1
lambda_U238 = 1.55125e-10 # Decay constant of uranium-238 in years^-1
lambda_Th232 = 4.9475e-11 # Decay constant of thorium-232 in years^-1

U_238_lambda = -1.54e-10
U_235_lambda = -9.72e-10
Th_232_lambda = -4.940e-11

U_238_initial_mol_g = U_238_mol_g / np.exp(-U_238_lambda * 4.5e9)
U_235_initial_mol_g = U_235_mol_g / np.exp(-U_235_lambda * 4.5e9)
Th_232_initial_mol_g = Th_mol_metal / np.exp(-Th_232_lambda * 4.5e9)

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


