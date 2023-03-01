import numpy as np
import matplotlib.pyplot as plt

Pressures = np.linspace(0, 100, 5)
colors = ['red', 'blue', 'green', 'purple', 'orange']

Temp = np.linspace(1000, 7000, 100) #K
# reduced the number of elements in Temp to make the code faster and use less memory

a_i_Ni = 0.46  #unitless
b_i_Ni = 2700 # K
c_i_Ni = -61 #K/ GPa


a_i_Si = 1.3 #unitless
b_i_Si = -13500 #K/GPa


a_i_O = .6 #Unitless
b_i_O =-3800 #K
c_i_O =22 #K/GPa

log_K_O_D_array = np.zeros((100, len(Pressures)))

for i, Pres in enumerate(Pressures):
    for j in range(100):
        log_K_O_D_array[j, i] = a_i_O + (b_i_O/Temp[j]) + (c_i_O*Pres)/(Temp[j])

plt.figure(1)
for i, Pres in enumerate(Pressures):
    plt.plot(1000/Temp, log_K_O_D_array[:, i], color=colors[i], label=f'log_K_OD at {Pres} GPa')

plt.xlabel("Temperature (K)")
plt.ylabel("log KOD")
plt.ylim(-5,2)
plt.xlim(min(1000/Temp),max(1000/Temp))
plt.legend()
plt.title("Logarithm of equilibrium constant for KOD")
plt.show()

plt.figure(2)
for i, Pres in enumerate(Pressures):
    log_K_Si_D = a_i_Si + (b_i_Si/Temp) 
    plt.plot(1000/Temp, log_K_Si_D, color=colors[i], label=f'log_K_SiD at {Pres} GPa')

plt.xlabel("Temperature (K)")
plt.ylabel("log KSiD")
plt.ylim(-5,2)
plt.xlim(min(1000/Temp),max(1000/Temp))
plt.legend()
plt.title("Logarithm of equilibrium constant for KSiD")
plt.show()

plt.figure(3)
for i, Pres in enumerate(Pressures):
    log_K_Ni_D = a_i_Ni + (b_i_Ni/Temp) + (c_i_Ni*Pres)/(Temp) 
    plt.plot(1000/Temp, log_K_Ni_D, color=colors[i], label=f'log_K_NiD at {Pres} GPa')

plt.xlabel("Temperature (K)")
plt.ylabel("log KNiD")
plt.ylim(-5,2)
plt.xlim(min(1000/Temp),max(1000/Temp))
plt.legend()
plt.title("Logarithm of equilibrium constant for KNiD")
plt.show()