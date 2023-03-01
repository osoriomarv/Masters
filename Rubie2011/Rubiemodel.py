# Rubie at al 2011 model
# Code will solve for the fo2 using a heteregeneous accretion model

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.integrate as integrate


#7 Guess an initial value of x'
#6 Determine y' from Eq. (S15)
#5 Obtain a' and b' from Eq. (S14a,b)
#4 Using Eq. (S17), determine z'
#3 Determine c' and d' from Eq. (S14c,d)
#2 Obtain  from Eq. (S6)
#1 Substitute determined values into Eq. (S11) and compare the resulting value of   with that determined from the model of Frost et al. (2010). Depending on the difference between the two values, adjust the value of x' and return to step 2, repeating the calculation (steps 2-7) iteratively until  obtained using Eq. (S11) converges with that determined from the model of Frost et al. (2010). Note that   has to be recalculated from the model of Frost et al. (2010) at each iteration.


# Define the composition

# wt%
Si = 10.68
Mg = 9.61
Fe = 18.43
Al = 1.04
Ca = 1.14
Ni = 1.077
Co = 0.0505
# Ppb
Nb = 300
Ta = 17.3
W = 110
# ppm
V = 60.3
Cr = 2646

# Partitioning of oxygen between metal Partitioning of oxygen between metal, magnesiowüstite and silicates 

#Define log Kd(si) = a - b/T + (c1*P + c2*P = c3P3)/T

def log_Kd(a,b,T,c1,c2,c3,P1,P2,P3):
    log_Kd = a - b/T + (c1*P1 + c2*P2 + c3*P3)/T
    return



# gamma_metal_FeO = activity of FeO 
# X_Metal_FeO = mole fraction of FeO in metal

#As described in the main text, the oxygen content of the bulk composition controls the oxygen fugacity. Oxygen partitioning is also critical for estimating the oxygen content of the Earth’s core. The thermodynamic treatment of oxygen partitioning between metal and silicate is derived from experimental data, obtained up to 70 GPa and 3500 K,  on the partitioning of oxygen between (i) liquid metal and magnesiowüstite, (ii) metal and silicate melt and (iii) the partitioning of FeO between magnesiowüstite and silicate melt (Asahara et al. 2007; Frost et al., 2010). The thermodynamic basis for the model considers the equilibrium between the FeO components in liquid metal and magnesiowüstite with the condition at equilibrium
u_metal_feo = U_mg_wustite_FeO
#where U_metal_feo is the chemical potentia of FeO in phase m. These chemical potentials are related to the mole fractions of FeO i each pashe, e.g. X_metal_FeO, through equations:
u_metal_feo = U_feo_P_T + R*T*X_metal_FeO + R*T*gamma_metal_FeO
#and
u_mg_wustite_feo = U_feo_P_T + R*T*X_mg_wustite_FeO + R*T*gamma_mg_wustite_FeO
#where U_mg_wustite_FeO is the standard state chemical potential, R is the gas constact and gamma_mg_wustite_FeO is, for the activity coefficient of FeO in the metallic phase. The standard state chemical potentials of FeO in liquid Fe and magnesiowüstite were determined using thermodynamic data extracted from the melting curve of FeO (Frost et al., 2010). The non-ideal mixing terms for FeO and Fe in the metallic liquid and magnesiowüstite are described using asymmetric Margules equations and are taken from fits to experimental data as described by Frost et al. (2010). 
#The partitioning of FeO between silicate liquid and liquid Fe-alloy is determined from the partitioning of FeO between magnesiowüstite and liquid Fe alloy using an empirical relation of Rubie et al. (2004) that was determined from the data of Trønnes (2000), Trønnes and Frost (2002) and Ito et al. (2004):

X_mg_wustite_FeO = 1.148*X_sil_FeO + 1.319*x_sil_fe0**2 

#Where x_sil_FeO is the mole fraction of FeO in the silicate liquid. this expression shows no dependance on P and T based on data obtained over 20 GPa pressure range


#S2. Solution of the mass balance to determine the compositions of co-existing silicate and metal liquids

# In this reaction:
# Fixed coefficients: 		
# x, y, z, a, b, c, d	(defined by starting compositions)
# Constant coefficients: 	
# u, m, n 		(defined by bulk composition)
# Unknown coefficients:	
# x_prime, y_prime, z_prime, a_prime, b_prime, c_prime, d_prime	(require 7 independent equations to solve).


# Four of the required seven equations are provided by the following mass balance equations
#Mass Balance for Fe, Ni, O and Si
#For Fe: x+a = x' + a'
#For Ni: y + b = y' + b'
#For O: x + y + 2z + c = x' + y' + 2z' + c'
#For Si: z + d = z' + d'

#mole_fraction
X_sil_FeO  = x_prime / (x_prime + y_prime + z_prime + (u + m + n))

#oxygen_partitioning 

K_D_O = (x_met_feo * x_met_o)/X_mg_wustite_FeO
#where
x_met_o = c_prime / (a_prime + b_prime + c_prime + d_prime)
#and
x_met_fe = a_prime / (a_prime + b_prime + c_prime + d_prime)

#and substituting these gives

K_O_D = (a_prime*c_prime)/(X_mg_wustite_FeO* (a_prime + b_prime + c_prime + d_prime)**2)

#For Ni Partitioning 

K_Ni_D = (X_met_ni*X_sil_fe)/(X_sil_ni*x_met_fe)

#so that

K_Ni_D = (((b_prime)/(a_prime+b_prime+x_prime+d_prime))*((x_prime/(x_prime+y_prime+z_prime+(u+m+n))))/((y_prime/(x_prime+y_prime+z_prime+(u+m+n)))*(a_prime/(a_prime+b_prime+c_prime+d_prime))))=x_prime*b_prime / y_prime*a_prime
#where K_ni_d is determined from a fit to the data of kegler et al

#for Si Partitioning 

#K_Si_D = X_met_si*X_sil_FeO**2/x_sil_sio2*X_met_fe**2

K_Si_D = ((d_prime/(a_prime + b_prime + c_prime + d_prime)*(x_prime/(x_prime+y_prime+z_prime+u+m+n))**2)/((z_prime/(x_prime+y_prime+z_prime+u+m+n))*(a_prime/(a_prime+b_prime+c_prime+d_prime))**2))
#Thus,
K_Si_D = x_prime**2 * d_prime*(a_prime+b_prime+c_prime+d_prime)/(a_prime**2 *z_prime*(x_prime+y_prime+z_prime+u+m+n))

K_Si_D = {see S1.3}


a_prime = (x+a-x_prime)
b_prime = (y+b - y_prime)
c_prime = (x + y + 2*z + c - x_prime - y_prime - 2*z_prime)
d_prime = (z + d - z_prime)

y_prime = x_prime*(y + b)/((x+a-x_prime)*K_Ni_D + x_prime)

K_Si_D = x_prime**2*(z+d-z_prime)*(a_prime+b_prime+x+y+3*z+c-x_prime-y_prime-3*z_prime+d)/(a_prime**2*z_prime*(x_prime+y_prime+z_prime+u+m+n))

#1:

#Set initial gues of x_prime

#Determine y_prime for T and P

y_prime = x_prime*(y + b)/((x+a-x_prime)*K_Ni_D + x_prime)

#Obtain a_prime + b_prime 

a_prime = (x + a - x_prime)
b_prime = (y + b - y_prime)

#Determine z_prime from the quadratic 
# [3*x_prime**2 - a_prime**2*K_Si_D]*z_prime**2 - [gamme*x_prime**2 + 3*alpha*x_prime**2 + a_prime**2*sigma*K_Si_D]*z_prime + alpha*gamma*x_prime**2 = 0

#Determine c_prime + d_prime 
c_prime = (x + y + 2*z + c - x_prime - y_prime - 2*z_prime)
d_prime = (z + d - z_prime)

# Obtain x_mg_wustite_FeO
X_mg_wustite_FeO = 1.148*X_sil_FeO + 1.319*X_sil_FeO**2

#7)	Substitute determined values into Eq. (S11) and compare the resulting value of   with that determined from the model of Frost et al. (2010).
#  Depending on the difference between the two values, adjust the value of x' and return to step 2, 
# repeating the calculation (steps 2-7) iteratively until  
# obtained using Eq. (S11) converges with that determined 
# from the model of Frost et al. (2010). 
# Note that   has to be recalculated from the model of Frost et al. (2010) at each iteration.

