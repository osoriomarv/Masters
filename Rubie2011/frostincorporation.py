import numpy as np
import pandas as pd
from Ironwustite import Calc_IW

import matplotlib.pyplot as plt

# # Oxidized Compostion
# x = 0.0672935926615133 # FeO
# y = 0 # NiO
# z = 0.2044532676017216 # SiO2
# u = 0.21257647049048262 # Mg
# m = 0.010312442368560392 # Al
# n = 0.015255276180658181 # Ca
# a = 0.4498580094305469 # Fe
# b = 0.04025094126651692 # Ni
# c = 0 # O
# d = 0 # Si

# x=  0.400966654
# y=	0.416897598
# z=	0.131973859
# u=	0.020224404
# m=	0.029918118
# n=	1.93658E-05
# a=  0.215330677
# b=	0.019266662
# c=	0.003302661
# d=	0

# #Reduced Compositon
# x = 0.7110226145582626
# y = 0.019134930618154617
# z = 1.5447274878339456
# u = 0.001704378306764529
# m = 0.5471971473544683
# n = 0.007702323063130593
# a = 0.05407015393419428
# b = 0.03154665921427567
# c = 0.08617081478221545
# d = 0.02018916906720182

#Reduced Composition
x=	0.000422782
y=	0
z=	0.411564281
u=	0.524869719
m=	0.025472964
n=	0.0376461
a=	0.422352054
b=	0.023479057
c=	0
d=	0.089968889

# #reduced composition
# x = 0.855
# y = 0.005
# z = 2.000
# a = 0.144
# b = 0.045
# c = 0.000
# d = 1.127
# u = 0.419
# m = 0.051
# n = 0.030

# Input parameters for a single temperature and pressure
Pressures = 10  # GPa

# Constants
tolerance = 1e-3  # set a tolerance for the convergence criterion
R = 8.314 # Gas Constant

# Create DataFrame for a_i, b_i, c_i values
data = {'Element': ['Ni', 'Si', 'O'],
        'a_i': [0.46, 1.3, 0.6],
        'b_i': [2700, -13500, -3800],
        'c_i': [-61, np.nan, 22]}  # use np.nan for missing value in Si

df = pd.DataFrame(data).set_index('Element').fillna(0)

# Extract a_i, b_i, c_i values from DataFrame using .loc[]
a_i, b_i, c_i = df.loc[['Ni', 'Si', 'O']].T.values

#Temp
Temp = np.linspace(2000,2001, 1)  # temperature range
max_iter = np.linspace(1,1000,1000)  # maximum number of iterations

# Initialize array for storing results
K_O_D_list = np.zeros(len(Temp))
K_O_D_frost_list= np.zeros(len(Temp))
X_mg_wustite_FeO_list = np.zeros(len(Temp))
del_iw_log_list = np.zeros(len(Temp))
fe_mols_list =np.zeros(len(Temp))
feo_mols_list=np.zeros(len(Temp))
x_prime_list = np.zeros(len(Temp))
a_prime_list = np.zeros(len(Temp))
Del_IW_list = np.zeros(len(Temp))
IW_FO2_list = np.zeros(len(Temp))
si_mols_met_list = np.zeros(len(Temp))

error_perc_list = np.zeros(len(max_iter))

ni_mols_met_list = np.zeros(len(Temp))
o_mols_met_list= np.zeros(len(Temp))

x_prime_iter_list = []

# Iterate over temperature range
x_prime = .4  # initial guess for x_prime at 1000 K and 10 GPa

for i in range(len(Temp)):

    log_K_O_D = a_i[2] + (b_i[2] / Temp[i]) + (c_i[2] * Pressures) / (Temp[i])

    log_K_Si_D = a_i[1] + (b_i[1] / Temp[i])

    log_K_Ni_D = a_i[0] + (b_i[0] / Temp[i]) + (c_i[0] * Pressures) / (Temp[i])

    # Convert from log to real units
    K_O_D_fischer = round(10 ** log_K_O_D, 5)
    for j in range(int(len(max_iter))):
        # Step 2
        y_prime = (x_prime * (y + b)) / ((x + a - x_prime) * (10 ** log_K_Ni_D) + x_prime)  # Eq S.15
        # Step 3
        a_prime = x + a - x_prime  # S.14 A
        print (a_prime + x_prime - a - x)
        b_prime = y + b - y_prime  # S.14 B
        print(b_prime + y_prime - b - y)
        # Step 4
        # Define alpha, gamma, sigma
        alpha = z + d
        gamma = a_prime + b_prime + x + y + 3 * z + c - x_prime - y_prime + d
        sigma = x_prime + y_prime + u + m + n
                # solve the quadratic equation to obtain z_prime
        A = 3 * x_prime ** 2 - a_prime ** 2 * (10 ** log_K_Si_D)  # S17 A
        B = gamma * x_prime ** 2 + 3 * alpha * x_prime ** 2 + a_prime ** 2 * sigma * (10 ** log_K_Si_D)  # S17 B
        C = alpha * gamma * x_prime ** 2  # S17 C
        z_prime = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A)  # S 17 Solve for z'

        # Step 5
        c_prime = x + y + 2 * z + c - x_prime - y_prime - 2 * z_prime  # S14 C
        print(c_prime + x_prime + y_prime + 2 * z_prime - c - x - y)
        #check mass balance for Fe
        d_prime = z + d - z_prime  # S14 D
        print(d_prime + z_prime - d - z)
        # Step 6
        X_sil_FeO = x_prime / (x_prime + y_prime + z_prime + (u + m + n))
        X_mg_wustite_FeO = (1.148 * X_sil_FeO )+ (1.319 * (X_sil_FeO ** 2))  # S6

        # Step 7
        # Sub into S11
        K_O_D_1 = a_prime * b_prime
        K_O_D_2 = X_mg_wustite_FeO * ((a_prime + b_prime + c_prime + d_prime) ** 2)
        K_O_D = K_O_D_1 / K_O_D_2
        
        #Incorporating Frost et al. 2010
        
        K_O_D_frost = (a_prime * c_prime) / (X_mg_wustite_FeO)
        
        # Check convergence
        # print(np.abs(K_O_D - K_O_D_fischer) / K_O_D_fischer)

        error_perc_list[j]=np.abs((K_O_D - K_O_D_frost) / K_O_D_frost)

        if abs((K_O_D - K_O_D_frost) / K_O_D_frost) < 0.1:
            #print(f"At Temp = {Temp[i]:.0f}K and P = {Pressures:.1f}GPa, K_O_D has converged to {K_O_D:.5f}")
            fo2 = 2*np.log10(X_mg_wustite_FeO/a_prime)
            feo_mols = x_prime/(x_prime + y_prime + z_prime + (u + m + n))
            fe_mols = a_prime/(a_prime + b_prime + c_prime + d_prime)
            ni_mols_met = b_prime/(a_prime + b_prime + c_prime + d_prime)
            o_mols_met = c_prime/(a_prime + b_prime + c_prime + d_prime)
            si_mols_met = d_prime/(a_prime + b_prime + c_prime + d_prime)
            del_iw_log = 2*np.log10(feo_mols/fe_mols)   
            IW_FO2 = Calc_IW(np.array([Temp[i]]))
            Del_IW = fo2 
            si_mols_met_list[i] = si_mols_met
            ni_mols_met_list[i] = ni_mols_met
            o_mols_met_list[i] = o_mols_met
            IW_FO2_list[i] = IW_FO2
            del_iw_log_list[i] = del_iw_log
            fe_mols_list[i] = fe_mols
            feo_mols_list[i] = feo_mols
            K_O_D_list[i] = K_O_D
            X_mg_wustite_FeO_list[i] = X_mg_wustite_FeO
            a_prime_list[i] = a_prime
            x_prime_list[i] = x_prime
            Del_IW_list[i] = Del_IW
            K_O_D_frost_list[i] = K_O_D_frost
            break
        else:
            x_prime = x - .1 #abs((K_O_D - K_O_D_frost) / K_O_D_frost) * x_prime + a_prime
            x_prime_iter_list.append(x_prime)
    else:
        #print(f"At Temp = {Temp[i]:.0f}K and P = {Pressures:.1f}GPa, K_O_D did not converge within 10 iterations.")
        K_O_D_list[i] = np.nan
        X_mg_wustite_FeO_list[i] = np.nan
        x_prime_list[i] = np.nan
        K_O_D_frost_list[i] = np.nan

plt.plot(range(1, len(error_perc_list) + 1), error_perc_list)
plt.xlabel('Iteration number')
plt.ylabel('Error percentage')
plt.title(f"Error percentage vs. iteration number at Temp = {Temp[i]:.0f}K and P = {Pressures:.1f}GPa")
plt.show()

# plt.plot(Temp, feo_mols_list, label='X_FeO')
# plt.plot(Temp, fe_mols_list, label='X_Fe')
# plt.legend()
# plt.xlabel('Temperature (K)')
# plt.ylabel('mol fraction')
# plt.text(1995, -2, 'Pressure: 10 GPa', fontsize=10) # add text annotation
# plt.show()

# # # create a list of NaN values to pad x_prime_iter_list
# # nan_list = np.repeat(np.nan, len(error_perc_list)-len(x_prime_iter_list))

# # # pad x_prime_iter_list with NaN values
# # x_prime_iter_list_padded = np.concatenate([x_prime_iter_list, nan_list])

# # # plot the two lists
# # plt.plot(max_iter, np.log10(error_perc_list))
# # plt.xlabel('Iteration')
# # plt.ylabel('error (%)')
# # plt.show()



# #print(f"K_O_D_fischer = {K_O_D_fischer:.5f}, K_O_D = {K_O_D:.5f}, X_mg_wustite_FeO = {X_mg_wustite_FeO:.5f}, x_prime = {x_prime:.5f}")

# # # Print the final results
# # print(f"K_O_D_list = {K_O_D_list}")
# # print(f"X_mg_wustite_FeO_list = {X_mg_wustite_FeO_list}")
# # print(f"a_prime = {a_prime_list}")
# # print(f"x_prime_list = {x_prime_list}")

# # # Plot Del_IW_list
# # plt.plot(Temp, del_iw_log_list, label='log fO2(IW) Reduced')
# # plt.legend()
# # plt.xlabel('Temperature (K)')
# # plt.ylabel('log fO2(IW)')
# # plt.text(1995, -2, 'Pressure: 10 GPa', fontsize=10) # add text annotation
# # plt.show()

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
# plt.plot(Temp, K_O_D_frost_list , label='K_O_D_frost')
# plt.legend()
# plt.xlabel('Temperature (K)')
# plt.ylabel('K_O_D')
# plt.show()
