import numpy as np

MW_FeO = 71.8464
MW_NiO = 74.6928
MW_SiO2 = 60.0843
MW_MgO = 40.3044
MW_Mg = 24.305
MW_Al2O3 = 101.9613
MW_Al = 26.9815386
MW_CaO = 56.0774
MW_Ca = 40.078
MW_Fe = 55.845
MW_Ni = 58.6934
MW_O = 15.9994
MW_Si = 28.0855

wt_percent = {'FeO': 0.06, 'NiO': 0.0, 'SiO2': 48.84, 'MgO': 41.79, 'Al2O3': 5.13, 'CaO': 4.17, 'Ni': 28.0/1000000, 'Fe': 85.59/100, 'O': 0.0, 'Si': 9.17}

mass_oxide = sum(wt_percent.values()) / 100.0  # total mass of oxides in grams
mass_silicate = 0.75  # grams of silicate
mass_metal = 0.25  # grams of metal
mass_element = {'Fe': 0.0006, 'Ni': 0.005, 'Si': 0.09175, 'Mg': 0.167325, 'Al': 0.038475, 'Ca': 0.031275, 'O': 0.388735}  # calculated from the given composition

for oxide in wt_percent:
    mass_element[oxide[:-1]] = mass_oxide * wt_percent[oxide] / eval('MW_'+oxide)  # convert oxides to elements

A = np.array([
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 1, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
    [1-MW_FeO/MW_Fe, 1-MW_NiO/MW_Ni, 1-MW_SiO2/MW_Si, 0, 0, 0, 1, 1, 1, 0], 
    [0, 0, 0, 1-MW_MgO/MW_Mg, 0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 1-MW_Al2O3/MW_Al, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 1-MW_CaO/MW_Ca, 0, 0, 0, 1], 
    [1, 0, 0, 0, 0, 0, -1+MW_FeO/MW_Fe, -MW_NiO/MW_Ni, -MW_SiO2/MW_Si, 0], 
    [0, 1, 0, 0, 0, 0, 0, -MW_NiO/MW_Ni, 0, 0], 
    [0, 0, 1, 0, 0, 0, 0, 0, -MW_SiO2/MW_Si, 0], 
    [0, 0, 0, 1, 0, 0, 0, 0, 0, -MW_MgO/MW_Mg], 
    [0, 0, 0, 0, 1, 0, 0, 0, 0, -MW_Al2O3/MW_Al], 
    [0, 0, 0, 0, 0, 1, 0, 0, 0, -MW_CaO/MW_Ca], 
    [-1, -1, -1, 0, 0, 0, 1-MW_FeO/MW_Fe, 1-MW_NiO/MW_Ni, 1-MW_SiO2/MW_Si, -MW_CaO/MW_Ca]
])

B = np.array([mass_element['Fe'], mass_element['Ni'], mass_element['Si'], mass_element['Mg'], mass_element['Al'], mass_element['Ca'], 0.75 - mass_silicate, 0, 0, 0, 0, 0, 0, 0, 0])

# solve for coefficients
coefficients = np.linalg.solve(A, B)

# print results
print('x =', coefficients[0])
print('y =', coefficients[1])
print('z =', coefficients[2])
print('u =', coefficients[3])
print('m =', coefficients[4])
print('n =', coefficients[5])
print('a =', coefficients[6])
print('b =', coefficients[7])
print('c =', coefficients[8])
print('d =', coefficients[9])



