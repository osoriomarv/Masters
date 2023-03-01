#Step 7 Build the frost model to see the difference between that which is calculted and that which
#Frost Predicts
import math 
def frost_KOD(T,P):

    # FeO Solid
    mu_standard = -27931+ 252.848*T - 46.12826*log(T) - .0057402984*(T**2)
    V_standard = 1.225
    alpha_TE = 3.481e-5 + 2.968e-9*T - 0.0806*(T**-2) - 0.0014437*(T**-1)
    K_t = 1500000 - 2.64e2*(T-298) + 0.01906*(T-298)**2
    K_prime = 4.305

    # FeO Liquid
    mu_standard_liq = mu_standard + 34008 - 20.969*T
    V_standard_liq = 1.3244
    alpha_TE_liq = 4.923e-5 + 2.968e-9*T - 0.0806*T**-2 - 0.001443*T**-1
    K_t_liq = 802655 - 100*(T-298)
    K_prime_liq = 4.397

return KOD
