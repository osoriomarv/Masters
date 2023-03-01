import pandas as pd
import numpy as np

def Calc_IW(T):
    # Load the data
    data = pd.read_excel('O2FeFeO.xlsx', sheet_name='Sheet3', usecols='B:D', skiprows=0, nrows=9).values

    A = data[0]
    B = data[1]
    C = data[2]
    D = data[3]
    E = data[4]
    F = data[5]
    G = data[6]
    H_Kj = data[7]
    S_298 = data[8]
    H_jmol = H_Kj * 1000


    # Temperature
    R = 8.3145
    T_ref = 298.15

    # Species
    num_T = len(T)
    num_a = len(A)

    # Preallocation
    t_inv = np.zeros(num_T)
    Cp_T = np.zeros((num_T, num_a))
    dH_species = np.zeros((num_T, num_a))
    S_Species = np.zeros((num_T, num_a))
    delHrxnlogFeO = np.zeros(num_T)
    delSrxnlogFeO = np.zeros(num_T)
    Grxn_FeO_Jmol = np.zeros(num_T)
    log_K_FeO = np.zeros(num_T)
    ln_K_FeO = np.zeros(num_T)
    IW_FO2 = np.zeros(num_T)
    IW_fo2 = np.zeros(num_T)
    Grxn_FeO_Jmol_IW = np.zeros(num_T)

    for i in range(num_T):
        # Inverse Temperature
        t_inv[i] = T[i] / 1000
        # Heat Capacity
        Cp_T[i, :] = A + B * t_inv[i] + C * t_inv[i]**2 + D * t_inv[i]**3 + E / t_inv[i]**2

        # Enthalpy and Entropy
        dH_species[i, :] = (A*t_inv[i]+ (B*t_inv[i]**2)/2 + (C*t_inv[i]**3)/3 + (D*t_inv[i]**4)/4 - E/t_inv[i]+ F - H_Kj)*1000
        S_Species[i, :] = (A*np.log(t_inv[i]) + B*t_inv[i] + C*t_inv[i]**2/2 + D*t_inv[i]**3/3 - E/(2*t_inv[i]**2) + G)

        # Fe + 0.5O2 = FeO
        Hrxn_29815_FeO = H_jmol[2] - (0.5 * H_jmol[1] + H_jmol[0])
        Srxn_29815_FeO = S_298[2] - (0.5 * S_298[1] + S_298[0])
        Grxn_29815_FeO = Hrxn_29815_FeO - 298.15 * Srxn_29815_FeO
        K_29815_feO = np.exp(-Grxn_29815_FeO / (R * T_ref))

        # Temperature and Pressure corrections
        Hrxn_FeO = Hrxn_29815_FeO + dH_species[i, 2] - (0.5 * dH_species[i, 1] + dH_species[i, 0])
        delHrxnlogFeO[i] = Hrxn_FeO  # Logs the Hrxn FeO calculations
        Srxn_FeO = S_Species[i, 2] - (0.5 * S_Species[i, 1] + S_Species[i, 0])
        delSrxnlogFeO[i] = Srxn_FeO  # Logs the entropy

        # IW buffer and Log_K Calc
        Grxn_FeO_Jmol_IW[i] = Hrxn_FeO - T[i] * Srxn_FeO
        log_K_FeO[i] = -Grxn_FeO_Jmol_IW[i] / (2.303 * R * T[i])
        ln_K_FeO[i] = -Grxn_FeO_Jmol_IW[i] / (R * T[i])
        IW_FO2[i] = np.exp(Grxn_FeO_Jmol_IW[i] / (-R * T[i])) ** -2
    return IW_FO2

