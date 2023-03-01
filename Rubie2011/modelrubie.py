import numpy as np
import math
import pandas as pd

#FeO = x, NiO = y, SiO2 = z, Mg = u, Al = m, Ca = n, Fe = a, Ni = b, Si = d
# U,M,M are kept the same 
#Calculate the coefficients for x, y, z, a, b, c, d
#FeO
x=0.00021306163760449946
#NiO
y=0
#SiO2
z = 0.20740848253616845
#Mg
u = 0.2645089404118509
#Al
m = 0.012837141286691067
#Ca
n = 0.018971812892312565
#Fe
a = 0.39100348811417207
#Ni
b = 0.021736352958884663
#O
c = 0
#Si
d = 0.08332072016231586

# Input parameters for a single temperature and pressure
Temp = 1000  # K
Pressures = 10  # GPa

# Constants
x_prime = 0.1628  # set an initial guess for x_prime
tolerance = 1e-5  # set a tolerance for the convergence criterion
max_iterations = 1  # set a maximum number of iterations
R = 8.314 # Gas Constant

# Create DataFrame for a_i, b_i, c_i values
data = {'Element': ['Ni', 'Si', 'O'],
        'a_i': [0.46, 1.3, 0.6],
        'b_i': [2700, -13500, -3800],
        'c_i': [-61, None, 22]}  # use None for missing value in Si

df = pd.DataFrame(data)
df = df.fillna(0)
df = df.set_index('Element')

# Extract a_i, b_i, c_i values from DataFrame using .loc[]
a_i_Ni, b_i_Ni, c_i_Ni = df.loc['Ni']
a_i_Si, b_i_Si, c_i_Si = df.loc['Si']
a_i_O, b_i_O, c_i_O = df.loc['O']

# Calculate the target difference
log_K_O_D = a_i_O + (b_i_O/Temp) + (c_i_O*Pressures)/(Temp)
log_K_Si_D = a_i_Si + (b_i_Si/Temp)
log_K_Ni_D = a_i_Ni + (b_i_Ni/Temp) + (c_i_Ni*Pressures)/(Temp)
# Convert from log to real units
K_O_D_fischer = 10**log_K_O_D
print(K_O_D_fischer)
K_Ni_D = 10**log_K_Ni_D
K_Si_D = 10**log_K_Si_D

target_difference = 0.05 * K_O_D_fischer
print(target_difference)

# Initialize arrays for K values
log_K_O_D = np.zeros(max_iterations)
log_K_Si_D = np.zeros(max_iterations)
log_K_Ni_D = np.zeros(max_iterations)

for i in range(max_iterations):

    # Step 2
    y_prime = x_prime * (y + b) / ((x + a - x_prime) * K_Ni_D + x_prime) #Eq S.15

    # Step 3
    a_prime = x + a - x_prime #S.14 A
    b_prime = y + b - y_prime #S.14 B

    # Step 4
    #Define alpha, gamma, sigma
    alpha = z + d
    gamma = a_prime + b_prime + x + y + 3*z + c - x_prime - y_prime + d
    sigma = x_prime + y_prime + u + m + n

    # solve the quadratic equation to obtain z_prime
    A = 3 * x_prime ** 2 - a_prime ** 2 * K_Si_D #S17 A
    B = gamma * x_prime ** 2 + 3 * alpha * x_prime ** 2 + a_prime ** 2 * sigma * K_Si_D #S17 B
    C = alpha * gamma * x_prime ** 2 #S17 C
    z_prime = (-B + np.sqrt(B ** 2 - 4 * A * C)) / (2 * A) # S 17 Solve for z'

    # Step 5
    c_prime = x + y + 2 * z + c - x_prime - y_prime - 2 * z_prime #S14 C
    d_prime = z + d - z_prime #S14 D

    # Step 6
    X_sil_FeO = x_prime / (x_prime + y_prime + z_prime + (u + m +n))
    X_mg_wustite_FeO = 1.148 * X_sil_FeO + 1.319 * X_sil_FeO ** 2 #S6 

    

    # Step 7
    #Sub into S11
    K_O_D_1 = a_prime * b_prime
    K_O_D_2 = X_mg_wustite_FeO * (a_prime + b_prime + c_prime + d_prime)**2

    K_O_D = K_O_D_1 / K_O_D_2 
    
    print(K_O_D)
    # Calculate the percentage difference between the model and target value
   
    difference = np.abs(np.log10(K_O_D) - 10**log_K_O_D)

    #print(percentage_difference)
    # Check if the percentage difference is below the tolerance
    if (difference < target_difference).all():
        break

    # Update x_prime based on the percentage difference
    x_prime = x_prime - 0.01 * (np.log10(K_O_D) - 10**log_K_O_D) / (10**log_K_O_D) * x_prime
    #x_prime = x_prime

    # print the values of difference and K_O_D at each iteration
    #print(f"Iteration {i}: difference={percentage_difference}, K_O_D={K_O_D}")


