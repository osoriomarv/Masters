import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#def NobleGasDistribution(SiO2,TiO2,Al2O3,FeO,Fe2O3,MgO,CaO,Na2O,K2O,radius,radius2,Taccretion):
    # Index the oxides to be used for calculations of molar fractions
#Funin = np.array([SiO2,TiO2,Al2O3,FeO,Fe2O3,MgO,CaO,Na2O,K2O])
Funin = np.array([48.88,1.73,16.77,9.71,0.00,6.65,9.86,3.62,1.93])


# R_v or radius vector. It takes in the initial radius and the radius2 to
# space create a spaced out vector in steps of 10.
#r_v = np.arange(radius, radius2+10, 10)
r_v = np.linspace(3310,7500,1000)


# Number of species considered
num_g = len(r_v)

# Preallocation
Volume = np.zeros(num_g)
r = np.zeros(num_g)
d = np.zeros(num_g)
den = np.zeros(num_g)
Pmass = np.zeros(num_g)
g = np.zeros(num_g)
SA = np.zeros(num_g)
Temp_bottom = np.zeros(num_g)
Pressure_bottom = np.zeros(num_g)
Matmhe = np.zeros(num_g)
mbse = np.zeros(num_g)
mcore = np.zeros(num_g)
px = np.zeros(num_g)
T_v = np.zeros(num_g)
p_v = np.zeros(num_g)

for j in range(num_g):
    
    # From Bouhifd et al. 2020
    # Dhe is the partitioning coefficient between Metal and Silicate
    Dhe = 0.017
    
    # Gravitational Constant - N*m^2/(kg^2);
    G = 6.67e-11
    
    # Kappa is refefered to as the grain oppacity from Mizuno 1980
    # Kappa is 1x10^-4 cm2/g but converted to m2/kg for unit purposes
    # Hayashi et al 1979 parameters are considered for pressure calculations
    Kappa = 1e-5 # m2/kg
    
    # Taacr or Time of Accretion is the amount of time by which the solar nebula
    # which encompasses a protoplanet mainly He-H (dominated by hydrogen)
    # Taccr=Taccretion;
    Taccr = 1e6
    
    # change data to r^cubed
    
    # Volume of planet Calculation in km^3
    Volume[j] = (4./3)*np.pi*(r_v[j]**3)
    
    # Density Calculation using a linear regression
    data = pd.read_excel('Planet_density.xlsx')
    radius_2 = data.iloc[:,1].values.astype(float)
    density = data.iloc[:,2].values.astype(float)

    
    # Polyfit will create a regression in which radius_2 meaning radius squared,
    # and density in g/cc will output a linear equation by which densities may
    # be acquired. radius will be squared in order to fit into the polyfit regression
    # which takes in the squared radius of four celestial bodies (Mars, Earth, Venus
    # and Vesta).
    p1 = np.polyfit(radius_2, density, 1)
    r[j] = r_v[j]**2
    
    # Regression using the radius squared conversion
    d[j] = p1[0]*r[j] + p1[1] # g/cm3
    
    # Unit conversion
    den[j] = d[j]*1
    # Calculation of Planetary Mass using volume and density 
    # Units in Kg
    Pmass[j] = den[j]*Volume[j] # kg
    
    # Calculation of gravity from calculated density and planetary mass
    # Gravitational Constant - N*m^2/(kg^2)
    # Mass - kg
    # radius in m^2
    # g in m/s^2
    g[j] = (G*Pmass[j])/(r_v[j]*1000)**2
    
    # ME is the mass of our current earth as reported from Sasaki and Nawazaka
    # 1990 - values can be obtained from wikipedia as well. We choose to report
    # the weight in Kg
    ME = 5.9722e24 # in Kg
    
    # Mean Molecular weight of gas from Hayashi et al 1979
    # From a hydrogen-helium cloud
    mu = 2.34
    
    # Temperature in Kelvin
    # Equation from Hayashi et al 1979. Temp_bottom is defined as the
    # temperature between a magma ocean and an atmosphere which has a
    # composition similar to a solar nebula
    Temp_bottom[j] = 4280*(mu/2.34)*(((den[j]/1e9)/4500)**(1/3))*(Pmass[j]/ME)**(2/3) # K
    
    # T_v will shorthand the Temp_bottom - a better definition would be to
    # create a temperature vector in kelvin
    T_v[j] = Temp_bottom[j]
    
    # Pressure_bottom = 2000; % bar
    # Pressure at the bottom calculation - calculation is under the same
    # conditions as temperature at the bottom. Pressure at the bottom is taken
    # from Sasaki and Nakazawa 1990 using the parameters provided by Mizuno et
    # al 1980 and Hayashi et al 1979.
    Pressure_bottom[j] = 9.02e7*((Kappa/1e-5)**-1)*(Taccr/1e6)*(mu/2.34)**4*((den[j]/1e9)/4500)*(Pmass[j]/ME)**3
    
    # Pressure at the bottom calculation. MPa 
    p_v[j] = Pressure_bottom[j]/1e6
    
    
    # From Elkin-Tanton 2008, a relation between partial pressure of the
    # surface, surface area and gravity of body. We assume that the partial
    # pressure of the He-H atmosphere is 20 percent of the whole atmosphere.
    px[j] = Pressure_bottom[j]*0.136 # pascal
    
    # Surface area of a spherical object using radius m^2 
    SA[j] = 4*np.pi*(r_v[j]*1000)**2
    
    # Elkins-Tanton 2008 shows an equation of to get a mass of a gas in the
    # atmosphere. Matm=surface area.*partial pressure/gravity of the body.
    
    # Calculation of mass of helium
    Matmhe[j] = ((SA[j]*px[j])/g[j])*1000 # grams
    
    # We assume that the Earth is 1/3 the core and 2/3 the mantle by weight.
    #
    # Our assumptions are based around the assumption that we are dealing with an
    # earthly body. Mass of the Bulk silicate earth in kg
    mbse[j] = (Pmass[j]*(2/3))*1000
    
    # Mass of core in kg
    mcore[j] = (Pmass[j]*(1/3))*1000
    
    # Calculation of Ionic Porosity from Iacono-Marziano et al. 2010;    
    # Return BSE
    BSE = Matmhe/1e6 # Convert Matmhe from grams to kg


plt.figure(10)
plt.plot(r_v,T_v,'-+',linewidth=1.5)
plt.xlabel('Radius Km', fontsize=14)
plt.ylabel('Temperature at bottom of atmosphere', fontsize=14)
plt.title('Temperature as radius increases')
plt.show()


import pandas as pd
import numpy as np

# Load data from Excel file
data4 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet1', usecols='D', skiprows=1).values
data3 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet4').values
data2 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet2', usecols='B', skiprows=1, nrows=10).values
data5 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet2', usecols='C', skiprows=1, nrows=6).values
data6 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet2', usecols='D', skiprows=1, nrows=6).values
data1 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet3').values
data = pd.read_excel('nbIP.xlsx', sheet_name='Sheet5').values

# Extract data columns
oxide_wt = data[:, 0]
oxide_8wt = data[:, 1]
V_1573 = data1[:, 0]
dvdt = data1[:, 1]
dvdt_1673 = data1[:, 2]
dPMVdP_dt = data1[:, 3]
Vs = data2[:, 0]
Lambda = data5[:, 0]
Kappa = data6[:, 0]
alpha = data3[:, 0]
beta = data3[:, 1]
Vca = data4[:, 0]

num_f = len(V_1573)
num_s = len(data2)

Temp_bottom = np.arange(1573, 1673 + 1)

pmv = np.zeros((num_f, len(Temp_bottom)))
rpmv = np.zeros((num_f, len(Temp_bottom)))
Molv_melt = np.zeros((num_f, len(Temp_bottom)))

# Loop over num_f
for x in range(num_f):
    Funin = data[x, 2:]
    Melt_total = np.sum(Funin)
    norm_element_comp = Funin * 100 / Melt_total
    mol_conversion = norm_element_comp / oxide_wt[x]
    molsum = np.sum(mol_conversion)
    mol_fraction = mol_conversion / molsum
    SiO2 = 0.25 * (mol_fraction[0] - 0.5 * (mol_fraction[3] + mol_fraction[5] + mol_fraction[6]) - mol_fraction[7] - mol_fraction[8])
    TiO2 = 0.25 * mol_fraction[1]
    Al2O3 = 0.375 * mol_fraction[2]
    FeO = 0.25 * mol_fraction[3]
    Fe2O3 = 0.375 * mol_fraction[4]
    MgO = 0.25 * mol_fraction[5]
    CaO = 0.25 * mol_fraction[6]
    Na2O = 0.375 * mol_fraction[7]
    K2O = 0.375 * mol_fraction[8]
    Funin2 = [SiO2, TiO2, Al2O3, FeO, Fe2O3, MgO, CaO, Na2O, K2O]
    mol8sum = np.sum(Funin2)
    mol8fraction = Funin2 / mol8sum
    
    for l in range(num_f):
        pmv[l] = dvdt_1673[l] + dPMVdP_dt[l] * (Temp_bottom - 1673)
        rpmv[l] = V_1573[l] + (dvdt[l] * (Temp_bottom - 1573)) + ((pmv[l] * 10 * (0.1 - 1)))
        wt_1 = mol_fraction * oxide_wt
        wt_onemol = np.sum(wt_1)
        wt8_1 = mol8fraction * oxide_8wt
        wt_8mol = np.sum(wt8_1)
        Molv_melt[l] = mol_fraction[l] * rpmv[l]
        sumMolv_melt = np.sum(Molv_melt)
        Mol8V_melt = sumMolv_melt * wt_8mol / wt_onemol
        
        for o in range(num_s):
            SiO21 = mol_fraction[0] * ((1 / (Temp_bottom - 273)) - (1 / 1300))
            Al2O31 = mol_fraction[2] * ((1 / (Temp_bottom - 273)) - (1 / 1300))
            MgOCaO = (mol_fraction[5] + mol_fraction[6]) * ((1 / (Temp_bottom - 273)) - (1 / 1300))
            Na2OK2O = (mol_fraction[7] + mol_fraction[8]) * ((1 / (Temp_bottom - 273)) - (1 / 1300))
            Fe2O31 = mol_fraction[4] * ((1 / (Temp_bottom - 273)) - (1 / 1300))
            Na2Oq = mol_fraction[7]**2 * ((1 / (Temp_bottom - 273)) - (1 / 1300))
            TSiO21 = SiO21 * Lambda[0]
            TAl2O31 = Al2O31 * Lambda[1]
            TMgOCaO = MgOCaO * Lambda[2]
            TNa2OK2O = Na2OK2O * Lambda[3]
            TFe2O31 = Fe2O31 * Lambda[4]
            TNa2Oq = Na2Oq * Lambda[5]
            TempLambda = TSiO21 + TAl2O31 + TMgOCaO + TNa2OK2O + TFe2O31 + TNa2Oq
            
            SiO22 = mol_fraction[0] * (Vs[o] * 10 + 0.1)
            Al2O32 = mol_fraction[2] * (Vs[o] * 10 + 0.1)
            MgOCaO2 = (mol_fraction[5] + mol_fraction[6]) * (Vs[o] * 10 + 0.1)
            Na2O1 = mol_fraction[7] * (Vs[o] * 10 + 0.1)
            FeOFe2O3 = (mol_fraction[3] + mol_fraction[4]) * (Vs[o] * 10 + 0.1)
            K2O1 = mol_fraction[8] * (Vs[o] * 10 + 0.1)
            PSiO22 = SiO22 * Kappa[0]
            PAl2O32 = Al2O32 * Kappa[1]
            PMgOCaO2 = MgOCaO2 * Kappa[2]
            PNa2O1 = Na2O1 * Kappa[3]
            PFeOFe2O3 = FeOFe2O3 * Kappa[4]
            PK2O1 = K2O1 * Kappa[5]

            Pressure = PSiO22 + PAl2O32 + PMgOCaO2 + PNa2O1 + PFeOFe2O3 + PK2O1
            
            for n in range(num_f):
                Vcamolf = mol_fraction[n] * Vca[n]
                sumVcamolf = np.sum(Vcamolf)
                Vsmolf = mol_fraction[n] * Vs[n]
                sumVsmolf = np.sum(Vsmolf)
            
                IPcoeff = (Mol8V_melt * wt_onemol) / wt_8mol
                IP = 100 - 100 / IPcoeff * (sumVcamolf + sumVsmolf + TempLambda + Pressure)
                
                negativelnxi = alpha[2] * IP + beta[2]


