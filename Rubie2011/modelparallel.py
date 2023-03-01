import numpy as np
import pandas as pd
from Ironwustite import Calc_IW
import matplotlib.pyplot as plt
import multiprocessing

#Reduced Composition
x=  0.000422782
y=  0
z=  0.411564281
u=  0.524869719
m=  0.025472964
n=  0.0376461
a=  0.422352054
b=  0.023479057
c=  0
d=  0.089968889

# Input parameters for a single pressure
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
temp_range = np.linspace(2000, 2002, 2)  # temperature range
max_iter = np.linspace(1,1000,1000)  # maximum number of iterations
# Initialize arrays for storing results
K_O_D_list = np.zeros(len(temp_range))
K_O_D_frost_list = np.zeros(len(temp_range))
X_mg_wustite_FeO_list = np.zeros(len(temp_range))
del_iw_log_list = np.zeros(len(temp_range))
fe_mols_list = np.zeros(len(temp_range))
feo_mols_list = np.zeros(len(temp_range))
x_prime_list = np.zeros(len(temp_range))
a_prime_list = np.zeros(len(temp_range))
Del_IW_list = np.zeros(len(temp_range))
IW_FO2_list = np.zeros(len(temp_range))
si_mols_met_list = np.zeros(len(temp_range))
ni_mols_met_list = np.zeros(len(temp_range))
o_mols_met_list = np.zeros(len(temp_range))
# Iterate over temperature range
x_prime = .4  # initial guess for x_prime at 1000 K and 10 GPa

for i in range(len(temp_range)):
    log_K_O_D = a_i[2] + (b_i[2] / temp_range[i]) + (c_i[2] * Pressures) / (temp_range[i])
    log_K_Si_D = a_i[1] + (b_i[1] / temp_range[i])
    log_K_Ni_D = a_i[0] + (b_i[0] / temp_range[i]) + (c_i[0] * Pressures) / (temp_range[i])

    # Convert from log to real units
    K_O_D_fischer = round(10 ** log_K_O_D, 5)

    # Initialize error percentage list
    error_perc_list = []

    for j in range(int(len(max_iter))):
        # Step 2
        y_prime = (x_prime * (y + b)) / ((x + a - x_prime) * (10 ** log_K_Ni_D) + x_prime)  # Eq S.15
        # Step 3
        # Step 3
        # Ensure that a_prime and b_prime are positive
        if x + a < x_prime:
            a_prime = 0
        else:
            a_prime = x + a - x_prime  # S.14 A
        if y < -b:
            b_prime = 0
        else:
            b_prime = y + b - y_prime  # S.14 B

        # Step 4
        # Define alpha, gamma, sigma
        alpha = z + d
        gamma = a_prime + b_prime + x + y + 3 * z + c - x_prime - y_prime + d
        sigma = x_prime + y_prime + u + m + n

            # solve the quadratic equation to obtain z_prime
        A = 3 * x_prime ** 2 - a_prime ** 2 * (10 ** log_K_Si_D)  # S17 A
        B = gamma * x_prime ** 2 + 3 * alpha * x_prime ** 2 + a_prime ** 2 * sigma * (10 ** log_K_Si_D)  # S17 B
        C = alpha * gamma * x_prime ** 2  # S17 C
        discriminant = B ** 2 - 4 * A * C
        if discriminant < 0:
            z_prime = 0
        else:
            z_prime = (-B + np.sqrt(discriminant)) / (2 * A)  # S 17 Solve for z'
        # Step 5
        c_prime = x + y + 2 * z + c - x_prime - y_prime - 2 * z_prime  # S14 C


        d_prime = z + d - z_prime  # S14 D
        sum_values = a_prime + b_prime + c_prime + d_prime + z_prime + x_prime + y_prime + u + m + n

        # Step 6
        X_sil_FeO = x_prime / (x_prime + y_prime + z_prime + (u + m + n))
        X_mg_wustite_FeO = 1.148 * X_sil_FeO + 1.319 * X_sil_FeO ** 2  # S6

        # Step 7
        # Sub into S11
        K_O_D_1 = a_prime * b_prime
        K_O_D_2 = X_mg_wustite_FeO * ((a_prime + b_prime + c_prime + d_prime) ** 2)
        K_O_D = K_O_D_1 / K_O_D_2

        # Incorporating Frost et al. 2010
        K_O_D_frost = (a_prime * c_prime) / (X_mg_wustite_FeO)

        # Append error percentage to list
        error_perc_list.append(np.abs((K_O_D - K_O_D_frost) / K_O_D_frost))
        print('x_prime = ', x_prime, 'y_prime = ', y_prime, 'z_prime = ', z_prime, 'a_prime = ', a_prime, 'b_prime = ',b_prime,'c_prime', c_prime, 'd_prime', d_prime, 'u = ', u, 'm = ', m, 'n = ', n)
        # Check convergence
        if abs((K_O_D - K_O_D_frost) / K_O_D_frost) < 0.01:
            # Convergence achieved
            print(f"K_O_D has converged to {K_O_D:.5f} and the value of K_O_D_frost is {K_O_D_frost:.5f}")
            print(f"At Temp = {temp_range[i]:.0f}K and P = {Pressures:.1f}GPa, K_O_D has converged to {K_O_D:.5f}")
            feo_mols = x_prime/(x_prime + y_prime + z_prime + (u + m + n))
            print(feo_mols)
            fe_mols = a_prime/(a_prime + b_prime + c_prime + d_prime)
            print(fe_mols)
            del_iw_log = 2 * np.log10(feo_mols / fe_mols)
            IW_FO2 = Calc_IW(np.array([temp_range[i]]))
            IW_FO2_list[i] = IW_FO2
            del_iw_log_list[i] = del_iw_log
            K_O_D_list[i] = K_O_D
            X_mg_wustite_FeO_list[i] = X_mg_wustite_FeO
            a_prime_list[i] = a_prime
            x_prime_list[i] = x_prime
            K_O_D_frost_list[i] = K_O_D_frost
            break
        elif j == int(len(max_iter)) - 1:
            # Maximum number of iterations reached without convergence
            print(f"At Temp = {temp_range[i]:.0f}K and P = {Pressures:.1f}GPa, K_O_D did not converge within {int(len(max_iter))} iterations.")
            K_O_D_list[i] = np.nan
            X_mg_wustite_FeO_list[i] = np.nan
            x_prime_list[i] = np.nan
            K_O_D_frost_list[i] = np.nan
            break
        else:
            # Update x_prime and continue iterating
            x_prime = x_prime - .01 * ((K_O_D - K_O_D_frost) / K_O_D_frost) * x_prime

        
# Plot error percentage vs. iteration number for the current temperature
plt.plot(range(1, len(error_perc_list) + 1), error_perc_list)
plt.xlabel('Iteration number')
plt.ylabel('Error percentage')
plt.title(f"Error percentage vs. iteration number at Temp = {temp_range[i]:.0f}K and P = {Pressures:.1f}GPa")
plt.show()

plt.plot(temp_range, fe_mols_list, label='fe_mols')
plt.plot(temp_range, feo_mols_list , label='feo_mols')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Mole fraction')
plt.show()