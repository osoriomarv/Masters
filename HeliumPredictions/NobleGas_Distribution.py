
def NobleGasDistribution(SiO2,TiO2,Al2O3,FeO,Fe2O3,MgO,CaO,Na2O,K2O,radius,radius2,Taccretion):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from math import log10, exp
        # Index the oxides to be used for calculations of molar fractions
    Funin = np.array([SiO2,TiO2,Al2O3,FeO,Fe2O3,MgO,CaO,Na2O,K2O])
    #Funin = np.array([48.88,1.73,16.77,9.71,0.00,6.65,9.86,3.62,1.93])


    # R_v or radius vector. It takes in the initial radius and the radius2 to
    # space create a spaced out vector in steps of 10.
    r_v = np.arange(radius, radius2+10, 10)
    #r_v = np.linspace(3310,7500,1000)

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
        Taccr=Taccretion
        #Taccr = 1e6
        
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


    # plt.figure(10)
    # plt.plot(r_v,T_v,'-+',linewidth=1.5)
    # plt.xlabel('Radius Km', fontsize=14)
    # plt.ylabel('Temperature at bottom of atmosphere', fontsize=14)
    # plt.title('Temperature as radius increases')
    # plt.show()


    # Read data from Excel sheets
    data4 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet1', usecols='D', nrows=10, engine='openpyxl').values
    data3 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet4', engine='openpyxl').values
    data2 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet2', usecols='B', nrows=10, engine='openpyxl').values
    data5 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet2', usecols='C', nrows=6, engine='openpyxl').values
    data6 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet2', usecols='D', nrows=6, engine='openpyxl').values
    data1 = pd.read_excel('nbIP.xlsx', sheet_name='Sheet3', engine='openpyxl').values
    data = pd.read_excel('nbIP.xlsx', sheet_name='Sheet5', engine='openpyxl').values

    # Variable definition
    oxide_wt = data[:, 1].astype(float)
    oxide_8wt = data[:, 2].astype(float)
    V_1573 = data1[:, 1].astype(float)
    dvdt = data1[:, 2].astype(float)
    dvdt_1673 = data1[:, 3].astype(float)
    dPMVdP_dt = data1[:, 4].astype(float)
    Vs = data2[:, 0].astype(float)
    Lambda = data5[:, 0].astype(float)
    Kappa = data6[:, 0].astype(float)
    alpha = data3[:, 1].astype(float)
    beta = data3[:, 2].astype(float)
    Vca = data4[:, 0].astype(float)

    # Number of species considered
    num_f = len(Funin)
    num_s = len(Lambda)
    num_a = len(alpha)

    # Preallocation
    norm_element_comp = np.zeros(num_f)
    mol_conversion = np.zeros(num_f)
    pmv = np.zeros((num_f, num_g))
    rpmv = np.zeros((num_f, num_g))
    wt_1 = np.zeros(num_f)
    wt8_1 = np.zeros(num_f)
    Molv_melt = np.zeros((num_f, num_g))
    Vcamolf = np.zeros(num_f)
    Vsmolf = np.zeros(num_f)

    # TempLambda = np.zeros(num_s)
    # PressureK = np.zeros(num_s)
    negativelnxi = np.zeros(num_a)

    # Preallocation for Oxide temperature and Pressure relations
    SiO21 = np.zeros(num_g)
    Al2O31 = np.zeros(num_g)
    MgOCaO = np.zeros(num_g)
    Na2OK2O = np.zeros(num_g)
    Fe2O31 = np.zeros(num_g)
    Na2Oq = np.zeros(num_g)
    TSiO21 = np.zeros(num_g)
    TAl2O31 = np.zeros(num_g)
    TMgOCaO = np.zeros(num_g)
    TNa2OK2O = np.zeros(num_g)
    TFe2O31 = np.zeros(num_g)
    TNa2Oq = np.zeros(num_g)
    TempLambda = np.zeros(num_g)

    # Pressure for oxides
    SiO22 = np.zeros(num_g)
    Al2O32 = np.zeros(num_g)
    MgOCaO2 = np.zeros(num_g)
    Na2O1 = np.zeros(num_g)
    FeOFe2O3 = np.zeros(num_g)
    K2O1 = np.zeros(num_g)
    PSiO22 = np.zeros(num_g)
    PAl2O32 = np.zeros(num_g)
    PMgOCaO2 = np.zeros(num_g)
    PNa2O1 = np.zeros(num_g)
    PFeOFe2O3 = np.zeros(num_g)
    PK2O1 = np.zeros(num_g)
    Pressure = np.zeros(num_g)

    # Preallocation of available Ionic Porosity
    IP = np.zeros(num_g)

    for x in range(num_f):
        Melt_total = sum(Funin)
        norm_element_comp[x] = (Funin[x]*100)/Melt_total
        mol_conversion[x] = norm_element_comp[x]/oxide_wt[x]
        molsum = np.sum(mol_conversion)
        mol_fraction = mol_conversion/molsum

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
        mol8fraction = Funin2/mol8sum

        for l in range(num_f):
            j = np.arange(len(Temp_bottom))

            pmv[l, j] = dvdt_1673[l] + dPMVdP_dt[l] * (T_v[j] - 1673)
            rpmv[l, j] = V_1573[l] + dvdt[l] * (T_v[j] - 1573) + pmv[l, j] * 10 * (p_v[j] - 0.1)

            wt_1[l] = mol_fraction[l] * oxide_wt[l]
            wt_onemol = np.sum(wt_1)
            wt8_1[l] = mol8fraction[l] * oxide_8wt[l]
            wt_8mol = np.sum(wt8_1)

            Molv_melt[l, j] = mol_fraction[l] * rpmv[l,j]

            sumMolv_melt = np.sum(Molv_melt)

            Mol8V_melt = sumMolv_melt * wt_8mol / wt_onemol

            for o in range(num_s):
                SiO21[j] = mol_fraction[0] * ((1 / (T_v[j] - 273)) - (1 / 1300))
                Al2O31[j] = mol_fraction[2] * ((1 / (T_v[j] - 273)) - (1 / 1300))
                MgOCaO[j] = (mol_fraction[5] + mol_fraction[6]) * ((1 / (T_v[j] - 273)) - (1 / 1300))
                Na2OK2O[j] = (mol_fraction[7] + mol_fraction[8]) * ((1 / (T_v[j] - 273)) - (1 / 1300))
                Fe2O31[j] = mol_fraction[4] * ((1 / (T_v[j] - 273)) - (1 / 1300))
                Na2Oq[j] = mol_fraction[7]**2 * ((1 / (T_v[j] - 273)) - (1 / 1300))
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
                    Vcamolf = np.zeros(num_f)
                    Vsmolf = np.zeros(num_f)
                    for m in range(num_f):
                        # Mol fraction multiplied by ionic molar volume calculated after
                        # prewitt 1969 assuming oxides are a sphere
                        # sumsVcamolf sums the individual oxides
                        Vcamolf[m] = mol_fraction[m] * Vca[m]
                        sumVcamolf = np.sum(Vcamolf)

                        # Vsmolf is the multiplication of mol fraction by Vs calculated from
                        # Iacono Marzianos Model
                        # SumVsmolf is the sum of all oxides
                        Vsmolf[m] = mol_fraction[m] * Vs[m]
                        sumVsmolf = np.sum(Vsmolf)

                    # Ip coeff is the molar volume times the weight based on one mol
                    # divided by the weight of melt based on 8 oxygens
                    IPcoeff = (Mol8V_melt * wt_onemol) / wt_8mol

                    # IP calculation based on IM equation which takes into account
                    # pressure and temperature calculations.
                    IP = 100 - 100 / IPcoeff * (sumVcamolf + sumVsmolf + TempLambda + Pressure)

                    # Negativelnxi is the natural log of the solubility constant.
                    # Alpha and beta are parameters derived from Carroll and Stolper with
                    # modification from IM.
                    negativelnxi = alpha[2] * IP + beta[2]


    tc = 5.1953
    pc_pa = 227460
    R = 8.314

    a = 0.42748 * (R ** 2 * tc ** 2.5) / pc_pa
    b = 0.08664 * (R * tc / pc_pa)

    RKEOS_g = (R * T_v[0]) / Pressure_bottom[0]

    num_search = 100
    RKEOS = np.zeros(num_search)

    V_mol = np.zeros(len(Pressure_bottom))
    A_squared = np.zeros(len(Pressure_bottom))
    B = np.zeros(len(Pressure_bottom))
    Z = np.zeros(len(Pressure_bottom))
    ln_fugacitycoeff = np.zeros(len(Pressure_bottom))
    fugacitycoeff = np.zeros(len(Pressure_bottom))
    ccSTPgHe_bar = np.zeros(len(Pressure_bottom))
    He_sol_g_g = np.zeros(len(px))
    solubility_he = np.zeros(len(Pressure_bottom))
    BseHe = np.zeros(len(Pressure_bottom))
    Cmetal = np.zeros(len(Pressure_bottom))
    core = np.zeros(len(Pressure_bottom))
    Hesum = np.zeros(len(Pressure_bottom))
    atmpercent = np.zeros(len(Pressure_bottom))
    BSE = np.zeros(len(Pressure_bottom))
    Core = np.zeros(len(Pressure_bottom))

    for p in range(len(Pressure_bottom)):
        v_v = np.logspace(log10(RKEOS_g / 10), log10(RKEOS_g * 2), num_search)
        for i in range(len(v_v)):
            RKEOS[i] = ((R * T_v[p]) / (v_v[i] - b)) - (a / ((v_v[i] * (v_v[i] + b)) * T_v[p] ** 0.5)) - Pressure_bottom[p]

        min_val_i = np.argmin(np.abs(RKEOS))

        V_mol[p] = v_v[min_val_i]
        RKEOS_g = V_mol[p]

        A_squared[p] = (0.42748 * tc ** (5 / 2)) / (pc_pa * T_v[p] ** (5 / 2))
        B[p] = 0.08664 * tc / (pc_pa * T_v[p])

        Z[p] = (Pressure_bottom[p] * V_mol[p]) / (R * T_v[p])

        ln_fugacitycoeff[p] = Z[p] - 1 - np.log(Z[p] - B[p] * Pressure_bottom[p]) - (A_squared[p] / B[p]) * np.log(1 + ((B[p] * Pressure_bottom[p]) / Z[p]))

        fugacitycoeff[p] = exp(ln_fugacitycoeff[p])

        ccSTPgHe_bar[p] = exp(negativelnxi[p]) * 22414 * (px[p] / 1e5) * (fugacitycoeff[p] / wt_onemol)

        He_sol_g_g[p] = (1e5 * (ccSTPgHe_bar[p] / 1e6)) / (298.15 * 8.314) * 4.0026

        BseHe[p] = He_sol_g_g[p] * (mbse[p])

        Cmetal[p] = Dhe * He_sol_g_g[p]

        core[p] = (Cmetal[p] * mcore[p])

        Hesum[p] = core[p] + BseHe[p] + Matmhe[p]

        for c in range(len(Hesum)):
            atmpercent[c] = (Matmhe[c] / Hesum[c]) * 100
            BSE[c] = (BseHe[c] / Hesum[c]) * 100
            Core[c] = (core[c] / Hesum[c]) * 100

    import matplotlib.pyplot as plt

    # Figure 1
    fig1, ax1 = plt.subplots()
    ax1.plot(r_v, atmpercent, color=[0, 0.5, 0], linewidth=1.5)
    ax1.plot(r_v, BSE, '-b', linewidth=1.5)
    ax1.plot(r_v, Core, '-k', linewidth=2.5)
    ax1.legend(['Atm Percent', 'Bulk Silicate Earth', 'Core'], loc='best', frameon=False)
    ax1.set_xlabel('Planetary Radius (Km)', fontsize=16)
    ax1.set_ylabel('% Helium Distribution', fontsize=16)
    ax1.set_title('Helium Distribution with Planetary Growth', fontsize=16)

    # Figure 2
    fig2, ax2 = plt.subplots()
    ax2.plot(r_v, core, '-k', linewidth=1.5)
    ax2.legend(['Core'], loc='best', frameon=False)
    ax2.set_xlabel('Planetary Radius (Km)', fontsize=16)
    ax2.set_ylabel('Log10 Helium Distribution (Grams)', fontsize=16)
    ax2.set_title('Helium Distribution with Planetary Growth', fontsize=16)
    ax2.set_yscale('log')

    # Figure 3
    fig3, ax3 = plt.subplots()
    ax3.plot(r_v, Matmhe, '-y', linewidth=1.5)
    ax3.plot(r_v, BseHe, '-k', linewidth=1.5)
    ax3.plot(r_v, Core, '-b', linewidth=1.5)
    ax3.legend(['Atm Percent', 'Bulk Silicate Earth', 'Core'], loc='best', frameon=False)
    ax3.set_xlabel('Planetary Radius (Km)', fontsize=16)
    ax3.set_ylabel('% Helium Distribution', fontsize=16)
    ax3.set_title('Helium Distribution with Planetary Growth', fontsize=16)

    # Figure 10
    fig10, ax10 = plt.subplots()
    ax10.plot(r_v, T_v, '-k', linewidth=1.5)
    ax10.set_xlabel('Radius Km', fontsize=14)
    ax10.set_ylabel('Temperature at bottom of atmosphere, K', fontsize=14)
    ax10.set_title('Temperature as radius increases', fontsize=16)

    # Figure 11
    fig11, ax11 = plt.subplots()
    ax11.plot(r_v, p_v, '-k', linewidth=1.5)
    ax11.set_xlabel('Radius, Km', fontsize=14)
    ax11.set_ylabel('Pressure at bottom of atmosphere, MPa', fontsize=14)
    ax11.set_title('Pressure as radius increases', fontsize=16)
    plt.show()

if __name__ == "__main__":
    NobleGasDistribution(48.88,1.73,16.77,9.71,0.00,6.65,9.86,3.62,1.93,3000,9000,1e6)

       


