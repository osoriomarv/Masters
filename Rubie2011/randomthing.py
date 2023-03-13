import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

start_time = time.time()

#Calculate Bulk composition defined by SiO2, MgO, FeO, Al2O3, CaO and NiO
#Based on Rubie et al 2011. Heterogenous accretion model.
# bulk composition,
# SiO2, MgO, FeO, Al2O3, CaO, NiO
melt_wt = [48.84, 41.79, 0.06, 5.13, 4.17, 0.0]
melt_mm = [28.19 + 16, 24.31 + 16, 55.85 + 16, 26.98 + 16 * 1.5, 40.08 + 16, 58.59 + 16]
melt_mols_i = np.array(melt_wt) / np.array(melt_mm)
melt_mol_frac_i = melt_mols_i / np.sum(melt_mols_i)
#Metal component
# Fe, Ni, O, Si
metal_wt = [85.89, 5.00, 0, 9.17]
metal_mm = [55.85, 58.59, 16, 28.19]
metal_mols_i = np.array(metal_wt) / np.array(metal_mm)
metal_mol_frac_i = metal_mols_i / np.sum(metal_mols_i)

#melt fraction needed to distribute metal vs silicate distribution
Fe_met_frac = 0.999
metal_silicate_ratio = Fe_met_frac * melt_mols_i[2] / (metal_mols_i[0] - Fe_met_frac * metal_mols_i[0])

metal_mols_i_scaled = metal_mols_i * metal_silicate_ratio

# bulk cation,Si,Mg,Fe,Al,Ca,Ni
bulk_cat = np.array([melt_mols_i[0] + metal_mols_i_scaled[3], melt_mols_i[1],
                     melt_mols_i[2] + metal_mols_i_scaled[0], melt_mols_i[3],
                     melt_mols_i[4], metal_mols_i_scaled[1]])
# bulk O
bulk_O = melt_mols_i[0] * 2 + melt_mols_i[1] + melt_mols_i[2] + melt_mols_i[3] * 1.5 + melt_mols_i[4]
bulk_comp = np.concatenate((bulk_cat, [bulk_O]))

#Silicate Component
x = melt_mols_i[2]
y = melt_mols_i[5]
z = melt_mols_i[0]
u = melt_mols_i[1]
m = melt_mols_i[2]
n = melt_mols_i[5]

#Metal component
a = metal_mols_i_scaled[0]
b = metal_mols_i_scaled[1]
c = metal_mols_i_scaled[2]
d = metal_mols_i_scaled[3]

# Input parameters for a single temperature and pressure
Pressures = 10  # GPa
Temp = np.linspace(2000,3500,1000)

# Search Parameters
tolerance = 1e-3  # set a tolerance for the convergence criterion
max_iter = np.linspace(1,1000,1000)  # maximum number of iterations

# Create DataFrame for a_i, b_i, c_i values
data = {'Element': ['Ni', 'Si', 'O'],
        'a_i': [0.46, 1.3, 0.6],
        'b_i': [2700, -13500, -3800],
        'c_i': [-61, np.nan, 22]}  # use np.nan for missing value in Si

df = pd.DataFrame(data).set_index('Element').fillna(0)

# Extract a_i, b_i, c_i values from DataFrame using .loc[]
a_i, b_i, c_i = df.loc[['Ni', 'Si', 'O']].T.values

# Initialize array for storing results
K_O_D_list = np.zeros(len(Temp))
X_mg_wustite_FeO_list = np.zeros(len(Temp))
del_iw_log_list = np.zeros(len(Temp))
fe_mols_list =np.zeros(len(Temp))
feo_mols_list=np.zeros(len(Temp))
x_prime_list = np.zeros(len(max_iter))
a_prime_list = np.zeros(len(Temp))
Del_IW_list = np.zeros(len(Temp))
IW_FO2_list = np.zeros(len(Temp))


x_prime_iter_list = []

# Iterate over temperature range
x_prime = .01  # initial guess for x_prime at 1000 K and 10 GPa

for i in range(len(Temp)):
    #Gas Constant
    R = 8.314 # Gas Constant
    #Oxygen
    log_K_O_D = (a_i[2] + (b_i[2] / Temp[i]) + (c_i[2] * Pressures) / (Temp[i]))
    K_O_D_fischer = round(10**log_K_O_D, 5)
    #Silicon
    log_K_Si_D = (a_i[1] + (b_i[1] / Temp[i]))
    K_Si_D = round(10**log_K_Si_D, 5)
    #Nickel
    log_K_Ni_D = (a_i[0] + (b_i[0] / Temp[i]) + (c_i[0] * Pressures) / (Temp[i]))
    K_Ni_D = round(10**log_K_Ni_D, 5)

    for j in range(int(len(max_iter))):
        # Step 2
        y_prime = (x_prime * (y + b)) / ((x + a - x_prime) * (K_Ni_D) + x_prime)  # Eq S.15
        # Step 3
        a_prime = x + a - x_prime  # S.14 A
        b_prime = y + b - y_prime  # S.14 B
        # Step 4
        # Define alpha, gamma, sigma
        alpha = z + d
        gamma = a_prime + b_prime + x + y + 3 * z + c - x_prime - y_prime + d
        sigma = x_prime + y_prime + u + m + n
        # solve the quadratic equation to obtain z_prime
        A = 3 * x_prime ** 2 - a_prime ** 2 * (K_Si_D)  # S17 A
        B = (gamma * x_prime ** 2 + 3 * alpha * x_prime ** 2 + a_prime ** 2 * sigma * (K_Si_D))  # S17 B
        C = alpha * gamma * x_prime ** 2  # S17 C
        z_prime = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A) # S 17 Solve for z'
        
        # Step 5
        c_prime = x + y + 2 * z + c - x_prime - y_prime - 2 * z_prime  # S14 C
        d_prime = z + d - z_prime  # S14 D

        #melt_mols = x_prime + y + z + u + m + n 

        melt_mols = x_prime + y_prime + z_prime  + u + m + n
        metal_mols = a_prime + b_prime + c_prime + d_prime

        print(melt_mols)
        # Step 6
        X_sil_FeO = x_prime / (x_prime + y_prime + z_prime + (u + m + n))
        X_mg_wustite_FeO = 1.148 * X_sil_FeO + 1.319 * X_sil_FeO ** 2  # S6

        FeO_mol_frac = x_prime / melt_mols
        NiO_mol_frac = y_prime / melt_mols
        SiO2_mol_frac = z_prime / melt_mols
        Fe_mol_frac = a_prime / melt_mols
        Ni_mol_frac = b_prime / melt_mols
        O_mol_frac = c_prime / melt_mols
        Si_mol_frac = d_prime / melt_mols

        melt_mols_list = melt_mols
        metal_mols_list = metal_mols

        # Step 7
        # Sub into S11
        K_O_D = Fe_mol_frac * O_mol_frac / FeO_mol_frac

        # Check convergence
        # print(np.abs(K_O_D - K_O_D_fischer) / K_O_D_fischer)

        if abs((K_O_D - K_O_D_fischer) / K_O_D_fischer).all() < 0.01:
            #print(f"At Temp = {Temp[i]:.0f}K and P = {Pressures:.1f}GPa, K_O_D has converged to {K_O_D:.5f}")
            del_iw_log = 2*np.log10(FeO_mol_frac/Fe_mol_frac)   
            
            K_O_D_list[i] = K_O_D
            X_mg_wustite_FeO_list[i] = X_mg_wustite_FeO
            a_prime_list[i] = a_prime
            x_prime_list[i] = x_prime
            break
        else:
            x_prime = x_prime + .01 * (K_O_D - K_O_D_fischer) / K_O_D_fischer * x_prime
            x_prime_iter_list.append(x_prime)
    else:
        #print(f"At Temp = {Temp[i]:.0f}K and P = {Pressures:.1f}GPa, K_O_D did not converge within 10 iterations.")
        K_O_D_list[i] = np.nan
        X_mg_wustite_FeO_list[i] = np.nan
        x_prime_list[i] = np.nan

end_time = time.time()
run_time = end_time - start_time

print(f"Total run time: {run_time} seconds")
# # Print the final results
# print(f"K_O_D_list = {K_O_D_list}")
# print(f"X_mg_wustite_FeO_list = {X_mg_wustite_FeO_list}")
# print(f"a_prime = {a_prime_list}")
# print(f"x_prime_list = {x_prime_list}")

# figure 1
fig1 = plt.figure(1)
plt.plot(max_iter, (K_O_D_fischer), ':k',)
plt.xlabel('FeO guess')
plt.ylabel('log10(KDO)')
plt.legend(['guess', 'Fischer'])
plt.show()

#figure 5
fig5 = plt.figure(5)
plt.plot(Temp, del_iw_log, '-k')
plt.xlabel('T, K')
plt.ylabel('IW')
plt.title('IW')
plt.show()



# # Plot Del_IW_list
# plt.plot(Temp, del_iw_log_list, label='log fO2(IW) Reduced')
# plt.legend()
# plt.xlabel('Temperature (K)')
# plt.ylabel('log fO2(IW)')
# plt.text(1995, -2, 'Pressure: 10 GPa', fontsize=10) # add text annotation
# plt.show()

# plt.plot(Temp, si_mols_met_list, label='Si')
# plt.plot(Temp, fe_mols_list, label='Fe')
# plt.plot(Temp, ni_mols_met_list, label='Ni')
# plt.plot(Temp, o_mols_met_list, label='O')
# plt.legend()
# plt.xlabel('Temperature (K)')
# plt.ylabel('Metal Mole Fraction')
# plt.text(1995, -2, 'Pressure: 10 GPa', fontsize=10) # add text annotation
# plt.show()

# # Plot K_O_D and K_O_D_fischer
# plt.plot(Temp, K_O_D_list, label='K_O_D')
# plt.plot(Temp, [round(10 ** (a_i[2] + (b_i[2] / T) + (c_i[2] * Pressures) / T), 5) for T in Temp], label='K_O_D_fischer')
# plt.legend()
# plt.xlabel('Temperature (K)')
# plt.ylabel('K_O_D')
# plt.show()








